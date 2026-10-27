# Findings

All checked issues from the review, in three groups by severity: S-001… (major), M-001… (minor),
N-001… (note). Within each group, issues are sorted by effort: quick, medium, design.
The format is described in `README.md`.

Next free IDs: S-049, M-074, N-030

## Summary

### Severe (major)

| ID | Severity | Confidence | Kind | Effort | Where | Title | Status |
|----|----------|------------|------|--------|-------|-------|--------|
| S-001 | major | confirmed | memory | quick | JacRes.cpp:1524 | `init_pres` without `fix_phase` reads `phRat[-1]` | fixed |
| S-002 | major | possible | physics | quick | constEq.cpp:587-592 | Softening length scale `Lm` looks inverted (Le/Lm instead of Lm/Le) | fixed (group C) |
| S-003 | major | confirmed | physics | quick | phase.cpp:1049 | Profile "Diabase-Caristan_1982" has activation volume Vn = 1 m^3/mol | fixed |
| S-004 | major | confirmed | logic | quick | surf.cpp:1444-1462 | Sediment model 2 uses values left over from the previous node | fixed |
| S-005 | major | confirmed | restart | quick | fastscape.cpp:2910, 2747 | FastScape restart loses basement and silt fraction if their output is off | open |
| S-006 | major | confirmed | logic | quick | fdstag.cpp:550-551 | Grid stretching moves the global box about 0, not about `bg_ref_point` | fixed |
| S-007 | major | confirmed | logic | quick | fdstag.cpp:727-737 | CFL time step ignores negative velocities on non-uniform grids | fixed |
| S-008 | major | confirmed | physics | quick | bc.cpp:723, 1087, 1121 | `bvel_temperature_constant` is scaled with the time unit, not as a temperature | open |
| S-009 | major | confirmed | logic | quick | bc.cpp:1680-1725, 2054 | User `bvel_velout` is always overwritten by the mass-balance value | fixed |
| S-010 | major | confirmed | logic | quick | bc.cpp:692-697, 1690-1700 | `velin_num_periods = 1` sets the inflow velocity to zero and ignores net periods | fixed |
| S-011 | major | confirmed | restart | quick | bc.cpp:1704-1715 | After restart, net-only inflow adds the net offset twice (function `static`) | fixed |
| S-012 | major | confirmed | physics | quick | bc.cpp:3125-3134 | 3D Gaussian plume inflow is not mass-conserving (missing radius² factor) | open |
| S-013 | major | confirmed | logic | quick | advect.cpp:1217-1228 | Garbage collection deletes the wrong markers when the AVD delete list is unsorted | fixed |
| S-014 | major | confirmed | parallel | quick | cvi.cpp:1305-1414 | STAGP interpolation treats processor boundaries as domain walls | fixed |
| S-015 | major | confirmed | logic | quick | marker.cpp:2358 | Oblique ridge axis is wrong unless the front edge is at y = 0 | open |
| S-016 | major | confirmed | input | quick | marker.cpp:984 | Layer `rand_ampl` is divided by the max phase ID instead of the length scale | open |
| S-017 | major | confirmed | input, memory | quick | marker.cpp:1996-2094 | Phase diagram load fails silently (density 0); size not checked against buffer | open |
| S-018 | major | confirmed | input, memory | quick | tools.cpp:627; marker.cpp:1687 | Polygon stretch factors after the last control polygon are left unset | fixed |
| S-019 | major | confirmed | logic | quick | advect.cpp:222-234; marker.cpp:210-260 | `rand_noiseGP` moves markers after cell mapping without remapping | fixed |
| S-020 | major | confirmed | logic | quick | phase_transition.cpp:756-763,811-812,870 | Box with `PhaseOutside = -1` (as documented) sets marker phase to -1 and stops the run | open |
| S-021 | major | confirmed | input | quick | phase_transition.cpp:85,100-119,227-259,667-711 | Typos in `Type`, `Parameter_transition`, `Name_Clapeyron` and bad `ID` are silently accepted | open |
| S-022 | major | confirmed | logic, memory | quick | dike.cpp:1312,1331 | Dike stress time averaging loops `i < sx+ny` instead of `sx+nx` | open |
| S-023 | major | confirmed | physics | quick | dike.cpp:463-479 | Dike heating scales the whole cell conductivity by the dike fraction | fixed |
| S-024 | major | confirmed | physics | quick | passive_tracer.cpp:560-561 | Passive tracer pressure and temperature are read half a cell too high | fixed |
| S-025 | major | confirmed | logic | quick | adjoint.cpp:2342-2470 | Observation points in the same cell overwrite each other's weights | open |
| S-026 | major | confirmed | logic | quick | adjoint.cpp:632-662 | Bounds: 0 means "no bound", and one bound pins every other parameter | open |
| S-027 | major | confirmed | memory | quick | adjoint.cpp:1159,1438,1534 | `fcconv` is written past its end (TAO not limited by `Inversion_maxit`; off by one without TAO) | open |
| S-028 | major | confirmed | memory | quick | objFunct.h:85; adjoint.cpp:1294-1450 | `factor2array[51]` is indexed up to 100 parameters | open |
| S-029 | major | confirmed | physics | medium | JacResAux.cpp:244-321 | Lithostatic and pore pressure are zero in boundary ghost cells, so boundary edges get half the value | fixed |
| S-030 | major | possible | physics | medium | JacResTemp.cpp:520-522 | Adiabatic heating sees a fake pressure jump at boundaries with inflow or outflow | open |
| S-031 | major | confirmed | logic | medium | constEq.cpp:1008-1060 | Phase diagram lookup reads one pressure row too low, and its edges are wrong | open |
| S-032 | major | confirmed | physics | medium | surf.cpp:1229-1253 | Slope-dependent erosion never erodes peaks and ridges (centered gradient) | open |
| S-033 | major | possible | physics | medium | surf.cpp:1203-1211 | Slope-erosion sub-step limit isn't a stability limit | open |
| S-034 | major | confirmed | logic | medium | LaMEMLib.cpp:765; fastscape.cpp:1599 | With FastScape, max-angle smoothing is undone every step | open |
| S-035 | major | confirmed | logic | medium | AVD.cpp:1192-1199 | AVD deletion ignores the Voronoi volumes (`tclaimed` is always 0) | fixed |
| S-036 | major | confirmed | physics | medium | phase_transition.cpp:592-594,667-711 | Clapeyron laws are wrong in SI units (C read as K, MPa/K read as Pa/K) | fixed (group C) |
| S-037 | major | confirmed | logic | medium | phase_transition.cpp:1303,1322,1334 | NotInAirBox uses only the first segment's `zbounds` | open |
| S-038 | major | confirmed | logic | medium | passive_tracer.cpp:585-607, 904-938, 1032-1049 | Passive tracers use a stale marker-to-cell list (after marker exchange) | open |
| S-039 | major | confirmed | input | medium | options.cpp:944-950, 126-162, 251-252 | Command-line options are ignored while solver defaults (MG levels, grid, thermal) are set | fixed |
| S-040 | major | confirmed | physics | medium | adjoint.cpp:2357-2380,1639-1663 | Point-velocity cost isn't the point misfit, and its adjoint gradient doesn't match it | open |
| S-041 | major | confirmed | physics | medium | adjoint.cpp:860-880,3455-3520 | Strain-rate targets aren't scaled, and the strain-rate derivative is incomplete | open |
| S-042 | major | confirmed | parallel | medium | adjoint.cpp:1744-1749 | `Adjoint_AdvectPoint` moves points on one CPU only, and again in every forward run | open |
| S-043 | major | possible | physics | design | fastscape.cpp:458-469, 646 | FastScape grid doesn't follow a stretching LaMEM box | open |
| S-044 | major | confirmed | parallel | design | bc.cpp:2595-2604, 2648-2657 | `fix_phase` / `fix_cell` lose faces on processor boundaries; results depend on CPU count | fixed (group C) |
| S-045 | major | confirmed | logic | design | advect.cpp:278,1086-1096 | With `exy_num_periods`, outflow markers are deleted, not wrapped | fixed |
| S-046 | major | possible | physics | design | LaMEMLib.cpp:715; tssolve.cpp:308 | Solve and marker advection use different dt when dt is not fixed | fixed (group C) |
| S-047 | major | confirmed | parallel, memory | design | adjoint.cpp:2505-2690,102-112 | Surface observations (`Ap = 3`) use wrong grid positions and crash with several CPUs in z | open |
| S-048 | major | possible | parallel | design | adjoint.cpp:1152-1157,1186-1189 | Parameter vectors have 100 entries per CPU, so TAO sees a vector that grows with CPU count | open |

### Minor

| ID | Severity | Confidence | Kind | Effort | Where | Title | Status |
|----|----------|------------|------|--------|-------|-------|--------|
| M-001 | minor | confirmed | logic | quick | JacResTemp.cpp:525,534 | Adiabatic heating efficiency is applied twice | open |
| M-002 | minor | confirmed | physics | quick | JacRes.cpp:558-602 | Top pressure shift isn't area-weighted and silently overwrites `p_shift` | open |
| M-003 | minor | confirmed | logic | quick | constEq.cpp:734-742 | `Phase_Melt` density ignores the `mfmax` cap | open |
| M-004 | minor | confirmed | logic | quick | phase.cpp:223-232 | Healing defaults are set only when printing is on | open |
| M-005 | minor | confirmed | logic | quick | surf.cpp:1531 | Sediment model 3 assumes the box starts at x = 0 | fixed |
| M-006 | minor | confirmed | input, parallel | quick | surf.cpp:1711-1754 | Initial surface perturbation: unscaled inputs, NaN with noise only | open |
| M-007 | minor | confirmed | physics | quick | surf.cpp:1651, 1662 | Topographic diffusion boundary is not zero-flux | open |
| M-008 | minor | confirmed | logic | quick | surf.cpp:1004, 1068, 1346 | Erosion and sedimentation levels overshoot by up to one step | fixed |
| M-009 | minor | possible | logic | quick | surf.cpp:461 | Surface outside the box silently gets zero velocity | fixed |
| M-010 | minor | confirmed | memory | quick | fastscape.cpp:1234-1317, 3100 | Memory leaks in FastScape code | open |
| M-011 | minor | possible | logic | quick | fastscape.cpp:1055 | Non-uniform FastScape interpolation can skip a node | open |
| M-012 | minor | confirmed | logic, memory | quick | tools.cpp:628 | `interpStretch` can leave stretch factors unset after the last control polygon | fixed |
| M-013 | minor | confirmed | input | quick | bc.cpp:848 | Phase 0 can't be used as permeable / plume mantle inflow phase | open |
| M-014 | minor | confirmed | logic | quick | bc.cpp:2097-2099, 2161-2163 | `relax_dist = 0` with `bvel_face_out` gives NaN when a cell center is exactly at top/bot | open |
| M-015 | minor | confirmed | physics | quick | bc.cpp:1480-1519, 3005-3017 | Plume temperature: different Gaussian width in BC and on injected markers; unset Tbot used as -1 | open |
| M-016 | minor | confirmed | input | quick | bc.cpp:675-679, 712-713, 1784-1802 | Typos in BC options are silently ignored; `exy` silently drops no-slip on vx | fixed |
| M-017 | minor | confirmed | logic | quick | matData.cpp:497-498, 529-530, 561-562 | Coarse edge viscosity at upper boundaries averages in zero ghost cells | fixed |
| M-018 | minor | confirmed | logic | quick | nlsolve.cpp:430-515 | Some SNES/KSP failure reasons are not printed or not treated as severe | fixed |
| M-019 | minor | possible | logic | quick | nlsolve.cpp:284 | MFFD mat-vecs count as SNES function calls; long Newton runs may hit max_funcs | won't fix |
| M-020 | minor | confirmed | memory | quick | advect.cpp:1485-1539 | ADVCheckCorners reads an empty buffer when the cell and its local neighbours are empty | fixed |
| M-021 | minor | confirmed | logic | quick | AVD.cpp:753-754,1136-1137 | AVD never refills empty volumes and at most doubles sparse ones | not an issue |
| M-022 | minor | confirmed | logic | quick | AVD.cpp:830-833 | `FindPointInCell` returns an error code as a cell index | fixed |
| M-023 | minor | confirmed | logic | quick | cvi.cpp:204-205 | RK2 deletes markers whose midpoint leaves the domain | fixed |
| M-024 | minor | confirmed | logic | quick | subgrid.cpp:673; LaMEMLib.cpp:731-734 | `ADVGetSedPhase` uses stale cell numbers when called for passive tracers | open |
| M-025 | minor | confirmed | memory | quick | marker.cpp:604-642 | Temperature-file setup loops over the wrong marker count; index can run out of range | fixed |
| M-026 | minor | confirmed | memory | quick | marker.cpp:294 | `ADVMarkSave` crashes when the save path has no '/'; APS not saved | partly fixed |
| M-027 | minor | confirmed | logic | quick | advect.cpp:1089 | Periodic x-shift uses `>` while the rank test uses `>=` | fixed |
| M-028 | minor | confirmed | physics | quick | advect.cpp:2053; bc.cpp:2952 | Adiabatic gradient uses \|z - Ztop\|, so rock above the reference level is heated | open |
| M-029 | minor | confirmed | input | quick | marker.cpp:1193-1200 | Ridge spreading rate can silently be 0 | open |
| M-030 | minor | confirmed | input | quick | marker.cpp:378; advect.cpp:156-180; cvi.cpp:1239 | Marker input not checked (negative phase, nmark_sub, nmin > nmax, avd <= 0, stagp_a = 0) | fixed |
| M-031 | minor | confirmed | memory | quick | cvi.cpp:800,826-828 | Memory leak of `vi->markind` in `ADVelReAllocStorage` | fixed |
| M-032 | minor | confirmed | memory | quick | dike.cpp:352,673 | Dikes without free surface read `phRat[-1]` | open |
| M-033 | minor | possible | logic | quick | dike.cpp:753-754,1095 | Dynamic dike can produce NaN (no neighbour within `filty`, or zero lithosphere thickness) | open |
| M-034 | minor | confirmed | restart | quick | dike.cpp:219,1588,1593 | Restart with dikes re-reads time stepping and resets the averaging counter | open |
| M-035 | minor | confirmed | input | quick | dike.cpp:212,306-345,1230 | Dike opening silently 0 without `bvel_velin`; `magPwidth` not scaled | open |
| M-036 | minor | confirmed | logic, memory | quick | passive_tracer.cpp:599, 370, 607, 1049 | Tracer melt phase search checks the wrong marker; empty search reads past the list | open |
| M-037 | minor | confirmed | input | quick | tssolve.cpp:56, 61 | CFL and CFLMAX range checks can never trigger | fixed |
| M-038 | minor | confirmed | parallel | quick | options.cpp:482 | Automatic coarse-grid CPU count uses the local, not the total, cell count | fixed |
| M-039 | minor | confirmed | logic | quick | paraViewOutPassiveTracers.cpp:32 | `out_ptr` overwrites the `Passive_Tracer` switch after tracers are set up | open |
| M-040 | minor | confirmed | memory | quick | paraViewOutBin.cpp:770-772 | `UpdatePVDFile` uses the file name after `free` in its error message | fixed |
| M-041 | minor | possible | restart | quick | paraViewOutBin.cpp:766-797 | After a restart from an older checkpoint, .pvd files keep old entries after `</VTKFile>` | fixed |
| M-042 | minor | confirmed | logic | quick | adjoint.cpp:920-947 | `Adjoint_ScaleCostFunction = Var` uses only the last point | open |
| M-043 | minor | confirmed | physics | quick | adjoint.cpp:3024-3040,3080-3095 | Scaling-law exponents are wrong for `log10` parameters | open |
| M-044 | minor | confirmed | logic | quick | adjoint.cpp:1866,1995,1328 | A parameter value of 0 gives a zero FD step (division by zero) | open |
| M-045 | minor | confirmed | parallel | medium | JacRes.cpp:1004-1017 | Dike strain-rate correction is missing in ghost cells, so results depend on CPU count | open |
| M-046 | minor | possible | physics | medium | constEq.cpp:798-822 | Dike opening still enters the stabilization stress and shear heating | fixed |
| M-047 | minor | possible | logic | medium | JacRes.cpp:1545-1703 | `init_lith_pres` only sets the pressure history, and its stop test is weak | fixed (group C) |
| M-048 | minor | confirmed | input | medium | phase.cpp:798-803; marker.cpp:523 | Phase temperature `T = 0` (°C) is treated as "not set" | not an issue |
| M-049 | minor | possible | physics | medium | constEq.cpp:603-626, 697, 756 | Shear and bulk moduli averaged with non-elastic phases counted as 0 | fixed (group C) |
| M-050 | minor | confirmed | logic | medium | surf.cpp:1167, 1223, 1643 | Slope erosion and topographic diffusion ignore periodic x | open |
| M-051 | minor | confirmed | parallel, restart | medium | tools.cpp:80, 293-322; fdstag.cpp:292-294 | Rank-0-only errors can hang; restart reads aren't checked | fixed |
| M-052 | minor | confirmed | input | medium | bc.cpp:2103-2123, 2167-2187 | Sign of `velin` on Right/Back faces flips depending on `bvel_face_out` | open |
| M-053 | minor | confirmed | parallel, input | medium | bc.cpp:1229-1274, 1148 | Fixed-cell file read: unchecked fread/stat, one-rank errors hang, file not closed on error | fixed |
| M-054 | minor | possible | input | medium | bc.cpp:2132-2142, 2213-2223 | CompensatingInflow has no mass-balance check (needs user `velbot`/`veltop`) | fixed |
| M-055 | minor | confirmed | logic | medium | matFree.cpp:637-684, 1304-1309, 1355-1360; matBFBT.cpp:208-209 | Matrix-free operator and diagonal treat the periodic x seam as a free-slip wall | fixed |
| M-056 | minor | confirmed | logic | medium | matData.cpp:642-645, 808-811 | Rediscretized coarse levels lose no-slip and pressure BCs on boundary ghosts | fixed |
| M-057 | minor | confirmed | parallel | medium | advect.cpp:1368-1568 | ADVCheckCorners depends on CPU count, reuses the same seed, misses markers on the centre plane | fixed |
| M-058 | minor | confirmed | parallel | medium | AVD.cpp:747-749,791-793 | AVD half-volume rule uses processor boundaries; nmax not halved | fixed |
| M-059 | minor | confirmed | logic | medium | AVD.cpp:813-816,865-884,944-946 | AVD "edge" control volumes are face volumes; printed labels swapped | fixed |
| M-060 | minor | confirmed | memory | medium | marker.cpp:1466,1622-1744 | Polygon file sizes not checked against buffers | fixed |
| M-061 | minor | confirmed | logic | medium | interpolate.cpp:102,158,220,377 | Corner output at the periodic x seam uses one side only | fixed |
| M-062 | minor | confirmed | physics | medium | marker.cpp:1013,1152,1215; phase_transition.cpp:1197,1339 | Halfspace temperatures use hard-coded kappa = 1e-6 (wrong in non-dim mode) | open |
| M-063 | minor | confirmed | parallel | medium | subgrid.cpp:85,606-608; advect.cpp:1906; tools.cpp:293,798 | Rank-local hard errors can hang; sedimentation-phase error on fast uplift | fixed |
| M-064 | minor | confirmed | logic | medium | phase_transition.cpp:883-905,922,1026-1029 | `ResetParam = APS` resets every step; side depends on type; broken for NotInAirBox | open |
| M-065 | minor | confirmed | logic | medium | phase_transition.cpp:851-864 | `BoxVicinity` uses y-size for z, and stale or zero bounds for NotInAirBox | open |
| M-066 | minor | confirmed | memory | medium | tssolve.cpp:381-423, 457-462 | Time-step schedule can overflow its 2000-entry arrays and cuts nstep_max | fixed |
| M-067 | minor | confirmed | logic | medium | LaMEMLib.cpp:713, phase_transition.cpp:967 | Time-step restart (CFLMAX) repeats phase transitions, so moving boxes move twice | open |
| M-068 | minor | confirmed | logic, input | medium | passive_tracer.cpp:262, 420-797, 63-91 | Passive tracer APS never updated, IDs wrong when nx != ny, typos silently mean "Always" | open |
| M-069 | minor | confirmed | physics | medium | outFunct.cpp:742-743, 767-768; JacResAux.cpp | SHmax/EHmax direction vectors are averaged to corners, which breaks near north-south | fixed |
| M-070 | minor | confirmed | memory | medium | adjoint.cpp:204-389 | `AllMaterialParameters` scan: writes past arrays, checks only the first block | open |
| M-071 | minor | confirmed | logic | medium | adjoint.cpp:1935,1947,1995-2024 | Adjoint solve leaves the forward KSP on the `as_` prefix and cells in a perturbed state | open |
| M-072 | minor | confirmed | cleanup | medium | adjoint.cpp:1123-1132,62-66 | `SyntheticForwardRun` saves a zero vector; the `OFdef = 0` path doesn't exist | open |
| M-073 | minor | confirmed | parallel | design | dike.cpp:111-120,959-1001,1314,1395 | Dynamic dike parallel limits are not checked up front; results depend on CPU layout | open |

### Notes

| ID | Severity | Confidence | Kind | Effort | Where | Title | Status |
|----|----------|------------|------|--------|-------|-------|--------|
| N-001 | note | confirmed | logic | quick | JacRes.cpp:917 | `dikeRHS` is not set when dikes are off | open |
| N-002 | note | confirmed | cleanup | quick | JacResAux.cpp:389; JacResTemp.cpp:85,107 | Material data changed inside compute loops; Controls copied for every cell | fixed |
| N-003 | note | confirmed | physics | quick | phase.cpp:1511; constEq.cpp:375,798 | `eta_min` is not a cutoff; it's added to every phase's viscosity | won't fix |
| N-004 | note | confirmed | input | quick | phase.cpp:526-548 | Creep input checks never fire; `Bd`/`Bn` are silently overwritten by `eta`/`eta0` | fixed (group C) |
| N-005 | note | confirmed | cleanup | quick | meltParam.cpp | Whole file is dead code (no callers) | open |
| N-006 | note | confirmed | logic | quick | phase.cpp:1471,1479 | PrintMatProp prints wrong values under some labels | open |
| N-007 | note | confirmed | cleanup | quick | surf.cpp, fastscape.cpp | Small clean-ups (duplicate reads, wrong comments/units, folder name rounding) | fixed |
| N-008 | note | confirmed | physics | quick | surf.cpp:288-292 | Average topography is a plain node average | fixed |
| N-009 | note | confirmed | cleanup | quick | Tensor.cpp:112-318, 493-603 | ISA / eigenvalue code is unused and has latent bugs | fixed |
| N-010 | note | confirmed | cleanup | quick | fdstag.cpp:80-83 | Non-unit bias is rejected, so all bias code is dead | fixed |
| N-011 | note | confirmed | input | quick | bc.cpp:2405-2411, 2234, 2369 | Velocity cylinder/box edge cases: NaN axis, skipped at initial guess unlike other BCs | fixed |
| N-012 | note | confirmed | cleanup | quick | bc.h, bc.cpp (various) | Wrong comments and messages, debug print every step, small clean-ups | open |
| N-013 | note | possible | logic | quick | AVD.cpp:1123-1125 | AVD point with no cells gives NaN coordinates | fixed |
| N-014 | note | confirmed | input | quick | input_file.dat:747-749; phase_transition.cpp:562-580 | Docs say Clapeyron parameters overwrite a predefined law; they are ignored | fixed |
| N-015 | note | confirmed | cleanup | quick | objFunct.cpp; adjoint.cpp:3684-4150 | Dead code: whole ObjFunct module and the pointwise FD kernels | fixed |
| N-016 | note | confirmed | physics | medium | JacResAux.cpp:357 | `gw_level_type = surf` gives a flat water table at the mean topography | not an issue |
| N-017 | note | possible | physics | medium | phase.cpp:888-1352 | Experiment type of the creep profiles not checked against the sources | not an issue |
| N-018 | note | confirmed | memory | medium | tools.cpp:608-609, 658-662 | Polygon helpers use stack arrays and can divide by zero | fixed |
| N-019 | note | confirmed | logic | medium | nlsolve.cpp:216-238 | Newton/Picard switch: divergence test uses initial residual; max-Newton gives one Picard step | fixed (group C) |
| N-020 | note | confirmed | physics | medium | lsolve.cpp, matrix.cpp, matAux.cpp | `jp_pgamma` means different things per solver; inv_eta Schur with penalty drops the eta term | not an issue |
| N-021 | note | confirmed | cleanup | medium | nlsolve, lsolve, multigrid, matData, matFree (various) | Option overrides, stale comments, dead test code, small clean-ups | partly fixed |
| N-022 | note | confirmed | cleanup | medium | marker, advect, AVD, cvi, subgrid (various) | Dead code, wrong comments and log lines, small clean-ups | fixed |
| N-023 | note | confirmed | cleanup | medium | phase_transition.cpp, dike.cpp (various) | Dead code, misleading names and comments, small clean-ups | open |
| N-024 | note | confirmed | input | medium | parsing.cpp (various) | Input parsing: silent bad numbers, octal integers, dropped values, rank-0 error hang | fixed |
| N-025 | note | confirmed | cleanup | medium | options.cpp, tssolve.cpp, scaling.cpp, LaMEMLib.cpp, docs | Output and time-step quirks, labels, docs | fixed |
| N-026 | note | confirmed | cleanup | medium | outFunct.cpp, paraViewOutBin.cpp, paraViewOutAVD.cpp, paraViewOutPassiveTracers.cpp | Output quirks: names, stale residual, dead vectors, unchecked inputs | fixed |
| N-027 | note | confirmed | input | medium | adjoint.cpp (several) | Adjoint input and small quirks | open |
| N-028 | note | confirmed | physics | design | JacRes.cpp:1086-1103 | FSSA only stabilizes the diagonal terms | won't fix |
| N-029 | note | confirmed | physics | design | subgrid.cpp:528-606 | Rock/air phase change keeps old stress, APS and T | won't fix |

## Details

## Severe (major)

### S-001 [major, confirmed] `init_pres` without `fix_phase` reads `phRat[-1]`
- Tags: major, confirmed, memory, quick
- Where: src/JacRes.cpp:1524 (commit 6124aa7d)
- What: JacResInitPres tests `svCell[iter].phRat[fixPhase] != 1.0`. `fixPhase` defaults to
  -1 (bc.cpp:574). Nothing requires `fix_phase` when `init_pres` is on. The other user of
  `fixPhase` (bc.cpp:2580) does check for -1.
- Why it matters: with `init_pres = 1` and no `fix_phase`, the first cell reads before the start
  of `svBuff` (out of bounds). Every other cell reads the last phase ratio of the previous
  cell. If that value is 1.0, the cell's pressure is not initialized. Cells that are skipped
  keep whatever was in the work vector from `DMGetGlobalVector`, which isn't zeroed.
- Fix idea: `if(fixPhase == -1 || svCell[iter].phRat[fixPhase] != 1.0)`, and `VecZeroEntries(gp)`.
- Status: fixed (group A)

### S-002 [major, possible] Softening length scale `Lm` looks inverted
- Tags: major, possible, physics, quick
- Where: src/constEq.cpp:587-592
- What: with `Lm` set, the strain limits are scaled by `Le/Lm`. For fracture-energy
  regularization, the softening should finish at a fixed displacement, and the band strain
  is about displacement/Le. So the critical strain should scale as `Lm/Le`. The code makes
  bigger cells soften over more strain, which is the opposite, and makes mesh dependence
  worse.
- Why it matters: anyone using `Lm` to reduce mesh dependence gets the reverse effect.
  There's no documentation, so check the intended definition with the author.
- Status: fixed (group C): confirmed by Anton (finer meshes must finish softening at larger
  strain); now `Lm/Le`, scaling documented in input_file.dat. No test uses `Lm`.

### S-003 [major, confirmed] Profile "Diabase-Caristan_1982" has Vn = 1 m^3/mol
- Tags: major, confirmed, physics, quick
- Where: src/phase.cpp:1049
- What: activation volume `m->Vn = 1;`, which is m^3/mol. Normal values are around 1e-5.
  At 1 GPa, P*V/RT is about 1e5, so `exp(-Q)` is 0 and dislocation creep is switched off
  everywhere except at the surface. All the other profiles use 0 or ~1e-5.
- Fix idea: probably `0` (the source in the comment gives no V). Check the reference.
- Status: fixed (group A)

### S-004 [major, confirmed] Sediment model 2 uses values left over from the previous node
- Tags: major, confirmed, logic, quick
- Where: src/surf.cpp:1444-1462
- What: `zprop` and `zpropn` are set only in some branches. If `dn<0`, only `zprop` is set and
  `zpropn` keeps the previous node's value. If `d<0`, it's the other way round. Then
  `dz = zpropn - zprop` mixes two nodes. At time 0, `z = zprop` can also be stale.
- Why it matters: wrong sediment thickness next to the margin, and the result depends on loop
  order, so also on how the domain is split across CPUs.
- Fix idea: set both values in every branch (probably both = InitLevel + hUp in the first two).
- Status: fixed (group A)

### S-005 [major, confirmed] FastScape restart loses basement and silt fraction if their output is off
- Tags: major, confirmed, restart, quick
- Where: src/fastscape.cpp:2910-2913 (copy), 2747-2755 (restart init), 3022/2977 (restart I/O)
- What: the basement and silt-fraction arrays are copied from FastScape only when
  `out_surf_basement` / `out_surf_silt_fraction` are on. But they are always written to the
  restart file and fed back with `fastscape_set_basement_` / `fastscape_init_f_`. With output
  off they are all zeros.
- Why it matters: after a restart, FastScape treats the whole topography as sediment (basement
  = 0), so kfsed/kdsed apply everywhere. The silt fraction is reset as well.
- Fix idea: always copy basement (and silt fraction if marine), or at least before writing a restart.
- Status: open

### S-006 [major, confirmed] Grid stretching moves the global box about 0, not about `bg_ref_point`
- Tags: major, confirmed, logic, quick
- Where: src/fdstag.cpp:550-551 (commit 6124aa7d)
- What: Discret1DStretch moves the nodes with x += eps*(x - ref), but the global bounds with
  `gcrdbeg *= (1+eps)` and `gcrdend *= (1+eps)`, which means ref = 0. When `bg_ref_point` isn't
  0 in a stretched direction, the box from FDSTAGGetGlobalBox moves away from the real grid by
  eps*ref every step, and the error adds up.
- Why it matters: the global box is used for surface clamping and erosion (surf.cpp), the
  top of the domain (bc.cpp:2949, advect.cpp:2044, phase_transition.cpp:1212, 1353),
  passive tracers, the FastScape ranges and the uniform-grid CFL step. No test hits this:
  t18 sets a ref point but uses simple shear only.
- Fix idea: gcrdbeg += eps*(gcrdbeg - ref); same for gcrdend.
- Status: fixed (group A); test with non-zero `bg_ref_point` not yet added

### S-007 [major, confirmed] CFL time step ignores negative velocities on non-uniform grids
- Tags: major, confirmed, logic, quick
- Where: src/fdstag.cpp:727-737 (called from advect.cpp:2006-2008)
- What: on non-uniform grids, idt = v/h is negative when v < 0, so it never raises idtmax.
  Only flow in the + direction limits the step. The uniform branch uses |v|. Any grid with
  segments of different cell size counts as non-uniform.
- Why it matters: when the fastest flow goes in the - direction (e.g. a sinking slab or
  sphere, where vz < 0), the time step can be far too large and markers can jump more than a
  cell. The CFLMAX restart check uses the same function, so it doesn't catch this either.
- Fix idea: idt = PetscAbsScalar(v)/h.
- Status: fixed (group A)

### S-008 [major, confirmed] `bvel_temperature_constant` is scaled with the time unit, not as a temperature
- Tags: major, confirmed, physics, quick
- Where: src/bc.cpp:723 (read), 1087 (print), 1121 (nondimensionalize) (commit 6124aa7d)
- What: the value is read with `scal->time` as scale, so it is divided by the time unit. The
  print multiplies it back, so the log shows the right number. Then line 1121 treats the
  divided value as °C: `(T/time_unit + Tshift)/T_unit`. Other temperatures use scale 1.0.
- Why it matters: with `bvel_temperature_inflow = Constant_T_inflow` in GEO or SI units the
  inflow markers get a wrong temperature (close to 0 °C in GEO units).
- Fix idea: read with scale 1.0 and print without `*scal->time`.
- Status: open

### S-009 [major, confirmed] User `bvel_velout` is always overwritten by the mass-balance value
- Tags: major, confirmed, logic, quick
- Where: src/bc.cpp:1680, 1699, 1719, 1725 (BCGetVelins), 2054 (BCApplyBoundVel)
- What: BCApplyBoundVel calls BCGetVelins first, and every branch of BCGetVelins sets
  `velout = -velin*(top-bot)/(bot-bz)`. The `bvel_velout` input read in BCCreate is never used.
- Why it matters: a documented input has no effect, silently.
- Fix idea: keep a flag for "user velout given" and only recompute when it is not set.
- Status: fixed (group B, batch 1)

### S-010 [major, confirmed] `velin_num_periods = 1` sets the inflow velocity to zero and ignores net periods
- Tags: major, confirmed, logic, quick
- Where: src/bc.cpp:692-697, 1690-1700
- What: `velin_array` is only read when `VelNumPeriods > 1`. BCGetVelins takes the branch
  `else if(bc->VelNumPeriods)` for any non-zero value, so with 1 it sets
  `velin = velin_array[0] = 0`. That branch also wins over the net-only branch, so
  `velin_net` periods are ignored too.
- Why it matters: a natural input ("1 period") gives zero inflow with no message.
- Fix idea: use `VelNumPeriods > 1` in that branch (and in the first test), or read the array
  for any value >= 1.
- Status: fixed (group A)

### S-011 [major, confirmed] After restart, net-only inflow adds the net offset twice
- Tags: major, confirmed, restart, quick
- Where: src/bc.cpp:1704-1715
- What: with only `velin_net` periods, the base inflow is stored in a function-`static`
  the first time BCGetVelins runs, then `bc->velin = base + net[kk]` is written back. The
  whole BCCtx is saved in the restart file, so after restart `bc->velin` already holds
  `base + net`, and the fresh `static` takes that as the base. Result: `base + 2*net`.
- Why it matters: restarted runs use a different inflow velocity than the original run.
- Fix idea: store the base value in BCCtx (set once in BCCreate) instead of a `static`.
- Status: fixed (group A)

### S-012 [major, confirmed] 3D Gaussian plume inflow is not mass-conserving (missing radius² factor)
- Tags: major, confirmed, physics, quick
- Where: src/bc.cpp:3125-3134
- What: the integral of `exp(-r²/c²)` over the box is `(π c²/4)·Δerf_x·Δerf_y`. The code uses
  `π/4·Δerf_x·Δerf_y/Area`, without `c²`. The 2D branch has the `c` factor and is right.
- Why it matters: with `Plume_VelocityType = Gaussian` and `Plume_Dimension = 3D` the outflow
  velocity is wrong, so net volume flux through the bottom is not zero. In an
  incompressible closed box this has no solution or forces flow through other boundaries.
- Fix idea: multiply `a, b, d, e` by `c*c`.
- Status: open

### S-013 [major, confirmed] Garbage collection deletes the wrong markers when the AVD delete list is unsorted
- Tags: major, confirmed, logic, quick
- Where: advect.cpp:1217-1228 (`ADVCollectGarbage`) (commit 6124aa7d)
- What: the loop `if(idel[ndel-1] != nummark-1) markers[idel[ndel-1]] = markers[nummark-1]; nummark--`
  is only right when `idel` is ascending. The exchange and cvi build it ascending, but
  `AVDDeletePointsMV` (AVD.cpp:1205) fills it by volume and Voronoi area.
  Example: nummark = 10, idel = [9, 3]: marker 9 is copied into slot 3 and survives,
  marker 8 is lost.
- Why it matters: when a rank deletes more markers in AVD control than it adds, markers that
  should go survive and unrelated markers vanish. Phases are silently lost or duplicated.
- Fix idea: sort `idel` before the loop, as `ADVCollectGarbageVec` in subgrid.cpp:419-420 now does.
- Status: fixed (group A)

### S-014 [major, confirmed] STAGP interpolation treats processor boundaries as domain walls
- Tags: major, confirmed, parallel, quick
- Where: cvi.cpp:1305/1314 (`IN == 0`, `IN == nx`), 1355/1364 (`JN`), 1405/1414 (`KN`) (commit 6124aa7d)
- What: `nx = fs->dsx.ncels` and `IN` are local. At an internal processor boundary the code
  uses the wall extrapolation `2*v[IN]-v[IN+1]` instead of the available ghost values.
- Why it matters: `interp = stagp` gives different marker paths for different CPU counts.
- Fix idea: test global indices (`sx+IN == 0`, `sx+IN == global nx`). The comment "ONLY FOR
  FREE SLIP BC" (cvi.cpp:1298) is not enforced either.
- Status: fixed 2026-10-26, b7a5c143 (cell-center velocities are computed first and exchanged with one
  ghost layer; also fixes a 2-cell interval at upper walls in the old code)

### S-015 [major, confirmed] Oblique ridge axis is wrong unless the front edge is at y = 0
- Tags: major, confirmed, logic, quick
- Where: marker.cpp:2358 (commit 6124aa7d)
- What: `x_oblique = (xL-xR)/(yF-yB) * y + xL;` should use `(y - yF)`. Divides by zero if
  yF == yB.
- Why it matters: plate age and halfspace temperature are wrong in the whole ridge box.
  Tests only use straight ridges.
- Fix idea: use `(y - y_ridgeFront)` and check yF != yB.
- Status: open

### S-016 [major, confirmed] Layer `rand_ampl` is divided by the max phase ID
- Tags: major, confirmed, input, quick
- Where: marker.cpp:984, 975, 2207 (commit 6124aa7d)
- What: `getScalarParam(..., "rand_ampl", ..., 1, (PetscScalar) maxPhaseID)` divides by the
  last argument, which should be `chLen`. With one phase (maxPhaseID = 0) it gives inf.
  `cosine` is range-checked against maxPhaseID too, so `cosine = 2` is silently ignored.
  The noise uses `rand()` per marker, so it depends on CPU count.
- Why it matters: layer interface noise has the wrong size, or inf/NaN layer bounds.
- Fix idea: pass `chLen`; check `cosine` against 1; use PetscRandom.
- Status: open

### S-017 [major, confirmed] Phase diagram load fails silently; size not checked against buffer
- Tags: major, confirmed, input, memory, quick
- Where: marker.cpp:1996-2000, 2090-2094, 2060-2088 (commit 6124aa7d)
- What: more than 8 diagrams ("buffer too small") or an unknown column count only print a
  message and return. The name is not stored, so constEq.cpp:990 sets `pd->rho = 0`.
  `nT*nP` is never checked against `_max_pd_sz_` (40100), so a 201×201 file writes past
  `rho_v`. `fopen`/`fscanf` results are not checked; printed max T/P are one step too wide.
- Why it matters: runs continue with zero density, or memory is overwritten.
- Fix idea: SETERRQ in both cases; check n <= `_max_pd_sz_`; check file I/O.
- Status: open

### S-018 [major, confirmed] Polygon stretch factors after the last control polygon are left unset
- Tags: major, confirmed, input, memory, quick
- Where: tools.cpp:627 (`interpStretch`); marker.cpp:1687-1696 (commit 6124aa7d)
- What: the tail is filled only if `CtrlPoly[numCtrlPoly-1] != numCtrlPoly-1`; it should
  compare with `numPoly-1`. E.g. one control polygon at position 1: all later slices use
  uninitialised stack values (only clamped to >= 0.1). In marker.cpp, PolyPos < 1, wrong
  order and unknown VolID are not checked; PolyID is read but unused.
- Why it matters: polygon bodies get random sizes.
- Fix idea: always fill up to numPoly-1; check 1 <= Pos <= Vol.num, increasing, VolID valid.
- Status: fixed (group A)

### S-019 [major, confirmed] `rand_noiseGP` moves markers after cell mapping without remapping
- Tags: major, confirmed, logic, quick
- Where: advect.cpp:222-234; marker.cpp:210-260 (commit 6124aa7d)
- What: `ADVMarkPerturb` runs after `ADVMapMarkToCells` and shifts markers by ±0.5 sub-cell.
  With `rand_noise` also on, markers can cross into the next cell or rank.
- Why it matters: `cellnum` is stale, so the first projection puts phase and T in the wrong
  cell; CheckMarkers misses it; markers off-rank cause "non-local marker" aborts.
- Fix idea: call `ADVMapMarkToCells` after the perturbation, or clamp to the cell.
- Status: fixed (group A): perturbed coordinate kept in host cell

### S-020 [major, confirmed] Box with `PhaseOutside = -1` sets marker phase to -1
- Tags: major, confirmed, logic, quick
- Where: src/phase_transition.cpp:756-763, 811-812, 870 (commit 6124aa7d)
- What: `nphc` is meant to be 0 when a transition does not change phases on both sides. But
  Box types never read `PhaseAbove/PhaseBelow`, which stay 0 (zeroed struct), so
  `PhaseAbove[0] >= 0 && PhaseBelow[0] >= 0` is always true and `nphc` is always 1. With
  `PhaseInside = 7`, `PhaseOutside = -1` (the documented way to "not change the phase
  outside", info/options/input_file.dat:760-762) and `BothWays` or `InsideToOutside`, a
  phase-7 marker outside the box gets `PH2 = PhaseOutside[0] = -1` and `P->phase = -1`.
  With `BoxVicinity = 1` only markers near the box are hit, otherwise all of them.
- Why it matters: `ADVCheckMarkPhases` then stops the run with "Detected markers with wrong
  phase". The documented example setup does exactly this as soon as phase-7 material leaves
  the box. No test uses `PhaseOutside = -1`.
- Fix idea: compute `nphc` from the arrays that belong to the type, and never assign a
  negative PH1/PH2 (keep `P->phase` instead).
- Status: open

### S-021 [major, confirmed] Typos in phase transition input are silently accepted
- Tags: major, confirmed, input, quick
- Where: src/phase_transition.cpp:85, 100-119, 227-259, 293, 667-711, 1411-1421
  (commit 6124aa7d)
- What:
  - Unknown `Type` (e.g. `Clapyron`): no error; `Type` stays 0 = `_Constant_` with
    `Parameter_transition = T` and `ConstantValue = 0`.
  - Unknown `Parameter_transition`: stays 0 = `_T_`, so the value is treated as a
    temperature; the check at line 293 passes.
  - Unknown `Name_Clapeyron`: `SetClapeyron_Eq` sets nothing, `neq = 0`, `Pres = {0,0}`,
    and the `neq != 1` branch gives `P >= 0` → every marker with positive pressure turns into
    the "above" phase.
  - `ID` uses `maxval = numPhtr-1`; with one transition `maxval = 0`, which disables the
    check, so `ID = 3` is accepted and the transition is silently never applied (only
    `matPhtr[0..numPhtr-1]` are looped). Same for `PhaseTransLinkLeft/Right` and dike `ID`
    (dike.cpp:163).
- Why it matters: a one-letter typo gives a model that runs but is wrong.
- Fix idea: add `else SETERRQ(...)` to each string chain; check `ID` range explicitly.
- Status: open

### S-022 [major, confirmed] Dike stress time averaging loops over the wrong range
- Tags: major, confirmed, logic, memory, quick
- Where: src/dike.cpp:1312, 1331 (commit 6124aa7d)
- What: after `DMDAGetCorners(DA_CELL_2D_tave, &sx, &sy, &sisc, &nx, &ny, ...)` the x loop
  is `for (i=sx; i < sx+ny; i++)`. It uses the local y count instead of nx.
- Why it matters: with `istep_nave > 1` (default 2) only the first `ny` columns are
  time-averaged; the rest keep the raw value, so dike relocation (`Set_dike_zones`) sees a
  mix. If `ny > nx` it writes outside the local array. The t26 dynamic-dike test uses
  `istep_nave = 1`, so this path is not tested.
- Fix idea: use `nx`.
- Status: open

### S-023 [major, confirmed] Dike heating scales the whole cell conductivity by the dike fraction
- Tags: major, confirmed, physics, quick
- Where: src/dike.cpp:463-479; JacResTemp.cpp:100,125 (commit 6124aa7d)
- What: `k` is the phase-averaged conductivity of the cell. `kfac` sums only
  `phRat[dike] * factor`, then `k = kfac*k`. A cell with 30 % dike and 70 % host rock gets
  `0.3*k` even below the solidus, where the factor is 1. With several dikes, `k` is
  multiplied again inside the loop with the running `kfac`.
- Why it matters: conductivity is far too low in every partly-dike cell (dike edges are
  always partly filled), which changes the thermal structure around the dike.
- Fix idea: `kfac = (1 - phRat[dike]) + phRat[dike]*factor`, applied once after the loop.
- Status: fixed (2026-10-23, Anton): as in the fix idea; jumps at T_sol/T_liq also smoothed

### S-024 [major, confirmed] Passive tracer pressure and temperature are read half a cell too high
- Tags: major, confirmed, physics, quick
- Where: src/passive_tracer.cpp:560-561 (commit 6124aa7d)
- What: P and T live in cell centers, but they are interpolated with the vz pattern:
  `InterpLin3D(p, II, JJ, K, ..., ccx, ccy, ncz)`. In z it uses the cell index `K` and node
  coordinates `ncz`, so the weight is computed against nodes while the values sit at centers.
  The result is the value at about `zp + dz/2`.
- Why it matters: P-T paths are the main use of passive tracers. With 1 km cells the pressure is
  off by about rho*g*dz/2 (~15 MPa) and the temperature by half a cell of the gradient (~10 K).
- Fix idea: use `KK` and `ccz` for p and T (same as for x and y).
- Status: fixed (group A)

### S-025 [major, confirmed] Observation points in the same cell overwrite each other's weights
- Tags: major, confirmed, logic, quick
- Where: src/adjoint.cpp:2342-2349, 2387-2394, 2433-2440 (and `xini` lines below each)
- What: projection weights are set with `=` in the local array. Two points that share nodes
  (same or neighbouring cell on one CPU) overwrite each other; only the last point counts there.
- Why it matters: with many points (e.g. dense surface data) part of the data is silently lost
  from the cost and the gradient.
- Fix idea: `+=` (the local-to-global step already uses ADD_VALUES), or one projection per point.
- Status: open

### S-026 [major, confirmed] Bounds: 0 means "no bound", and one bound pins every other parameter
- Tags: major, confirmed, logic, quick
- Where: src/adjoint.cpp:632-662; used at 1052-1055
- What: `if (ts)` treats a bound of 0 as "not given". A missing side is set to `InitGuess`.
  Any bound sets `Ab = 1`, which turns on TAO BLMVM with these vectors for all parameters. So:
  a parameter with no bounds gets `[InitGuess, InitGuess]` and can't move; one with only an
  upper bound can't go below its start value; `LowerBound = 0` is ignored. Parameters from
  `AllMaterialParameters` (index < iStart) never get bounds set, so they are `[0, 0]` and TAO
  pushes them to 0.
- Why it matters: inversions silently freeze parameters as soon as one bound is used.
- Fix idea: default bounds to `-PETSC_INFINITY` / `PETSC_INFINITY`; use the "found" flag, not
  the value, to detect a given bound.
- Status: open

### S-027 [major, confirmed] `fcconv` is written past its end
- Tags: major, confirmed, memory, quick
- Where: src/adjoint.cpp:1159 (size), 1438 (built-in descent), 1534 (TAO), 1095-1098 (print)
- What: `fcconv` has `Inversion_maxit` entries (default 50) per CPU. With TAO, `count` grows by
  one per function evaluation and is never limited by `Inversion_maxit` (no
  `TaoSetMaximumIterations`; the "1500" check only prints). Evaluation 51 writes past the end.
  Without TAO, `count` reaches `maxit` before the break, so `fcconvar[maxit]` is one past the end.
  The summary loop reads the same entries.
- Why it matters: heap corruption in normal inversion runs.
- Fix idea: pass `maxit` to TAO and size the array `maxit + 1`, or keep the history in a
  `std::vector`.
- Status: open

### S-028 [major, confirmed] `factor2array[51]` is indexed up to 100 parameters
- Tags: major, confirmed, memory, quick
- Where: src/objFunct.h:85; src/adjoint.cpp:1294, 1325, 1346, 1408, 1449
- What: the line-search factors are a fixed array of 51, but loops run to `mdN` (up to
  `_MAX_PAR_ = 100`).
- Why it matters: with more than 51 parameters and `Inversion_EmployTAO = 0`, writes past the
  array into the next `ModParam` fields (`maxfac`, `vel_scale`, ...).
- Fix idea: size it `_MAX_PAR_`.
- Status: open

### S-029 [major, confirmed] Lithostatic and pore pressure are zero in boundary ghost cells
- Tags: major, confirmed, physics, medium
- Where: src/JacResAux.cpp:244-321 (lithostatic), :339-411 (pore); used at
  src/JacRes.cpp:1195-1198, 1302-1305, 1410-1413
- What: `lp_lith` and `lp_pore` are zeroed and then filled only for owned cells.
  `LOCAL_TO_LOCAL` fills ghosts between processors but not the ghosts outside the domain.
  The edge loops in JacResFormResidual average 4 surrounding cells, including those
  outside-domain ghosts. Pressure and temperature don't have this problem: their ghosts are
  set by two-point constraints and FDSTAGSetEdgeCornerCenter.
- Why it matters: on edges at the domain boundary, `pc_lith` and `pc_pore` come out at about
  half the real value (a quarter at corners). Edge viscosity (activation volume, with
  `p_litho_visc` on by default) and yield stress (`p_litho_plast`, pore pressure) are then
  wrong at the boundary. This matters most for no-slip walls, where boundary edge stresses
  are not zero.
- Fix idea: after computing them, fill the outside ghosts by copying the nearest interior
  cell (a zero-gradient copy, same as SET_TPC with no constraint), then call
  FDSTAGSetEdgeCornerCenter.
- Status: fixed (group B, batch 3)

### S-030 [major, possible] Adiabatic heating sees a fake pressure jump at boundaries
- Tags: major, possible, physics, medium
- Where: src/JacResTemp.cpp:520-522
- What: the adiabatic term uses lithostatic pressure differences with the neighbour cell,
  including outside-domain ghosts, which are 0 (same cause as S-029). At a boundary cell,
  `(Pc - 0)/bdx * vx` is used as the gradient.
- Why it matters: wherever the velocity across the boundary is not zero (inflow or outflow
  side walls, bottom inflow, permeable sides), the cells next to the boundary get a large
  fake heating or cooling term ~ alpha*T*P*v/dx. With no flow through the walls the term is
  zero, so it's harmless there.
- Fix idea: the same ghost fill as S-029. Also, the one-sided boundary term should use
  zero gradient.
- Status: open (try t37 or a model with inflow boundaries, with `Adiabatic_Heat = 1`)

### S-031 [major, confirmed] Phase diagram lookup reads one pressure row too low; edges wrong
- Tags: major, confirmed, logic, medium
- Where: src/constEq.cpp:1008-1060 (setDataPhaseDiagram); loader at src/marker.cpp:2055-2088
- What: the table is stored with T changing fastest (checked in test/t09_PhaseDiagrams/*.in),
  so the row for (iP, iT) is `nT*iP + iT`, with iP and iT starting at 0. The code uses
  `nT*(indP-1) + indT`, so it reads values one pressure step (dP) lower than the weights
  assume. Also:
  - `if(indT[0] < 1)` and `if(indP[0] < 1)` also catch index 0, a valid one. Inside the
    first interval there's no interpolation; the code uses the value at index 1.
  - Upper clamps set `indT[1] = nT`, one past the last column, which wraps into the next
    P row. `ind[3] = nT*nP` is one past the end of the table. It doesn't crash only
    because the buffer is bigger (40100) and zero-filled, so density is pulled toward 0
    at the high-P/high-T corner.
  - `maxT/maxP = min + n*d` are printed one step too high (should be n-1).
- Why it matters: density and melt fraction from phase diagrams lag by one pressure step
  (about 90 MPa, roughly 3 km, for the test diagram) and are wrong near table edges. Regression
  tests were probably made with this code, so they won't notice.
- Fix idea: `ind = nT*indP + indT`, clamp the indices to [0, n-2] with weights clamped to
  [0, 1], then re-baseline t09, t16, t21 and t23.
- Status: open

### S-032 [major, confirmed] Slope-dependent erosion never erodes peaks and ridges
- Tags: major, confirmed, physics, medium
- Where: src/surf.cpp:1229-1253 (also the slope estimate at 1164-1190)
- What: the slope is a centered difference, (z[i+1]-z[i-1])/(2dx). At a local peak or on a
  ridge crest this is 0, so the node isn't eroded, while its flanks are. Example: z = [0, 1, 0]
  gives slope 0 in the middle. A zig-zag (odd-even) pattern is also invisible to it.
- Why it matters: peaks grow relative to their surroundings instead of wearing down, which is
  the opposite of what the model should do. New feature (PR #80), covered only by t38.
- Fix idea: use an upwind (Godunov or Osher) gradient for this Hamilton-Jacobi equation, e.g.
  |grad z|^2 = max(D-x,0)^2 + min(D+x,0)^2 + same for y.
- Status: open

### S-033 [major, possible] Slope-erosion sub-step limit isn't a stability limit
- Tags: major, possible, physics, medium
- Where: src/surf.cpp:1203-1211
- What: sub-steps keep the drop per step below 0.5*dx_min, i.e. k*S^n*dt <= 0.5*dx. For n = 1
  the front moves at speed k whatever the slope, so a stable explicit step needs k*dt/dx <= 1.
  The current rule allows k*dt/dx up to 0.5/S_max, which is above 1 for slopes below 0.5
  (most real topography).
- Why it matters: in gentle terrain the update can overshoot and oscillate. Together with S-032
  this can create spikes. Not reproduced yet.
- Fix idea: use dt <= C*dx_min/(n*k*S_max^(n-1)) (for n >= 1), together with an upwind gradient.
- Status: open

### S-034 [major, confirmed] With FastScape, max-angle smoothing is undone every step
- Tags: major, confirmed, logic, medium
- Where: src/LaMEMLib.cpp:765; src/fastscape.cpp:1599 (also 1617, 1651, 1722, ...)
- What: after FastScapeRun, LaMEM smooths `surf->gtopo` (surf_max_angle). On the next step,
  FastScapeRun starts again from its own saved topography (`gtopo_fs`, `gtopo_refine`, ...), not
  from `surf->gtopo`. So the smoothing is lost, and LaMEM and FastScape work with different
  surfaces.
- Why it matters: the smoothing does nothing lasting, and the phase ratios use a surface that
  FastScape never sees. t37 uses surf_max_angle = 45.
- Fix idea: smooth before handing the topography back, or copy the smoothed surface into the
  FastScape arrays (or disallow surf_max_angle with surf_mode = 2).
- Status: open

### S-035 [major, confirmed] AVD deletion ignores the Voronoi volumes
- Tags: major, confirmed, logic, medium
- Where: AVD.cpp:1192-1199 (commit 6124aa7d)
- What: `area[i] = A->chain[i].tclaimed;` but `tclaimed` is only counted in
  `AVDInjectPointsMV` (AVD.cpp:1110). In the delete path it is still 0 from `AVDCreate`, so
  the sort is a no-op and the first `n - nmax` markers are deleted.
- Why it matters: the whole AVD computation is wasted for deletion; markers in sparse parts of
  a crowded volume can be removed and leave holes.
- Fix idea: count claimed cells per point before sorting (share the code with the inject path).
- Status: fixed (AVD batch)

### S-036 [major, confirmed] Clapeyron laws are wrong in SI units
- Tags: major, confirmed, physics, medium
- Where: src/phase_transition.cpp:592-594, 667-711; scaling.cpp:157,172 (commit 6124aa7d)
- What: slopes are divided by `scal->stress/scal->temperature`. In GEO units `stress` is in
  MPa, so the input is MPa/K as documented. In SI units `stress` is in Pa, so the same number
  is read as Pa/K (10^6 too small). Also `T0` is converted with `Tshift`, which is 0 in SI,
  so the predefined `T0` values (in °C, e.g. 1537) are taken as Kelvin (273 K off).
- Why it matters: every predefined Clapeyron law (and user values following the [MPa/K]
  comment) is wrong in SI setups.
- Fix idea: store predefined laws in SI (Pa/K, K) and convert with `stress_si` and an
  explicit °C→K shift; document input units per unit system.
- Status: fixed (group C)

### S-037 [major, confirmed] NotInAirBox uses only the first segment's `zbounds`
- Tags: major, confirmed, logic, medium
- Where: src/phase_transition.cpp:1303, 1322-1323, 1334 (commit 6124aa7d)
- What: `zbounds` is read for every segment (`2*nsegs` values) but the inside test and the
  temperature profiles always use `zbounds[0]` and `zbounds[1]`.
- Why it matters: segments with a different depth range (as in the documented example,
  input_file.dat:781) get the wrong vertical extent.
- Fix idea: store per-cell z bounds like `celly_xboundL/R`, or document that only one z range
  is supported and read only two values.
- Status: open

### S-038 [major, confirmed] Passive tracers use a stale marker-to-cell list
- Tags: major, confirmed, logic, medium
- Where: src/passive_tracer.cpp:585-607, 1032-1049, 904-938; call order LaMEMLib.cpp:725-734
  (commit 6124aa7d)
- What: `ADVAdvectPassiveTracer` runs after `ADVAdvect` and `ADVExchange`, but before `ADVRemap`,
  which is where `ADVMapMarkToCells` rebuilds `markstart`/`markind`. `ADVCollectGarbage` moves
  markers into the holes, so the old indices point to other markers (or past `nummark`). If the
  marker storage grew, `ADVReAllocStorage` gives a fresh `markind` filled with zeros, so every
  cell points to marker 0. The cell IDs are also from the grid before `BCStretchGrid`.
- Why it matters: the "closest marker" phase used for melting tracers, delayed activation and
  numerical sedimentation of tracers comes from the wrong marker. Tracer phase output is wrong.
- Fix idea: call the tracer update after `ADVRemap`, or call `ADVMapMarkToCells` first.
- Status: open

### S-039 [major, confirmed] Command-line options are ignored while solver defaults are set
- Tags: major, confirmed, input, medium
- Where: src/options.cpp:944-950 (setSolverOptions), 126-162, 251-252 (commit 6124aa7d)
- What: `setSolverOptions` saves and clears all PETSc options, then runs
  `solverOptionsSetDefaults`, and only then puts the command-line options back. Inside, the code
  builds a scaling object and a temporary grid (`FDSTAGCreate`) and reads `act_temp_diff` and
  `act_steady_temp`. All of these read the input file only, so `-nel_x`, `-cpu_x`, `-nel_z`,
  `-act_temp_diff` etc. given on the command line are not seen.
- Why it matters: the number of multigrid levels and the block counts are computed for a different
  grid than the one that runs. The run may fail at MG setup (too many levels) or use poor
  settings. Thermal solver options (`ts_*`) are not set when temperature is turned on from the
  command line. Command-line overrides are common in parameter sweeps and in LaMEM.jl.
- Fix idea: read the solver-option inputs before clearing, or clear only after the defaults are
  computed (insert defaults first, then re-insert the saved options as now).
- Status: fixed (group B): no options clear in setSolverOptions

### S-040 [major, confirmed] Point-velocity cost isn't the point misfit, and its adjoint gradient doesn't match it
- Tags: major, confirmed, physics, medium
- Where: src/adjoint.cpp:2357-2380 (and Vy/Vz copies), 1639-1663 (commit 6124aa7d)
- What: for `Adjoint_ObservationPoints = 1` the target vector `xini` is built by scaling the
  current nodal velocities: `xini_k = target * v_k / v_point`. The cost is then
  `1/2 sum_k w_k (v_k - xini_k)^2` over the 8 nodes, and `dF/dx = w * (x - xini)` treats `xini`
  as fixed, though it depends on `x`. This only equals `1/2 (v_point - target)^2` when the 8
  nodal velocities are the same.
- Reproduced (debug build, t08 CompareGradients.dat with `CostFunction`, point 0.31 0.43 0.37):
  cost 8.80e-4 vs. true 8.68e-4; gradient eta[1]: FD -4.4e-10, adjoint +4.2e-9 (wrong sign);
  eta[0]: FD 1.72e-5, adjoint 1.49e-5. At the symmetric point 0.5 0.5 0.5 both agree.
- Also: a zero nodal or interpolated velocity gives inf/NaN in `xini` (division by `v_k`, `f_k`).
- Why it matters: inversions (TAO or the built-in descent) get a gradient that doesn't belong
  to the cost; line searches can fail or converge to the wrong place. No test covers this mode.
- Fix idea: `F = 1/2 sum_i (P_i . x - target_i)^2`, `dF/dx = sum_i P_i (P_i . x - target_i)`.
- Status: open

### S-041 [major, confirmed] Strain-rate targets aren't scaled, and the strain-rate derivative is incomplete
- Tags: major, confirmed, physics, medium
- Where: src/adjoint.cpp:860-880, 3455-3520; JacRes.cpp:661-677
- What:
  - `Value` for Exx..E2nd is stored as given (`Ae = ts`), but compared with non-dimensional
    strain rates (`Parameter - Ae`). The printed model value is converted to 1/s, so the user
    must give 1/s, which the cost then treats as non-dimensional.
  - The observed quantity is the effective deviatoric strain rate from JacResGetEffStrainRate,
    which includes the elastic history term `h*I2Gdt`, not the actual strain rate.
  - `dExx/dv` uses only `±1/dx` on Vx. Exx is deviatoric, so it is `2/3·(±1/dx)` on Vx and
    `-1/3·(±1/dy)`, `-1/3·(±1/dz)` on Vy, Vz. Same for Eyy, Ezz. E2nd has no derivative at all.
- Why it matters: strain-rate inversions minimise the wrong misfit; gradients are off.
- Fix idea: divide `Ae` by `scal.strain_rate`; use the total strain rate; add the missing terms
  (or forbid E2nd).
- Status: open

### S-042 [major, confirmed] `Adjoint_AdvectPoint` moves points on one CPU only, and again in every forward run
- Tags: major, confirmed, parallel, medium
- Where: src/adjoint.cpp:1744-1749
- What: `Ax/Ay/Az += v*dt` runs only on the CPU that holds the point, so the CPUs disagree on
  where it is; when it crosses a CPU border it can be lost or counted twice. The coordinates
  live in `IOparam`, which survives across `LaMEMLibMain` calls, so every FD run, adjoint run and
  inversion iteration moves the points further.
- Why it matters: with `Adv = 1`, results depend on CPU count and on how many forward runs were
  done; FD and adjoint gradients use different point positions.
- Fix idea: keep the start positions and reset them at the start of each forward run; reduce
  the moved positions over all CPUs.
- Status: open

### S-043 [major, possible] FastScape grid doesn't follow a stretching LaMEM box
- Tags: major, possible, physics, design
- Where: src/fastscape.cpp:458-469 (ranges set once), 646 (mode 0 "update range" never called)
- What: rangeX/Y/Z and the FastScape grid are set once at start-up. With background strain
  (BCStretchGrid), LaMEM's box changes, but topography is still mapped node by node and
  clamped to the old z range. The velocities sent to FastScape also include the stretching part.
- Why it matters: with bg strain + FastScape, topography ends up at the wrong x/y position and
  is clamped to wrong limits. Not tested; t37 has no background strain.
- Fix idea: update the range each step (mode 0), or reject bg strain with surf_mode = 2.
- Status: open

### S-044 [major, confirmed] `fix_phase` / `fix_cell` lose faces on processor boundaries
- Tags: major, confirmed, parallel, design
- Where: src/bc.cpp:2595-2604 (BCApplyPhase), 2648-2657 (BCApplyCells)
- What: for each owned cell the code also sets the face at `i+1`, `j+1`, `k+1`. For the last
  owned cell in a direction, that face belongs to the next rank and is only a ghost here.
  BCApply then calls LOCAL_TO_LOCAL (INSERT_VALUES), which overwrites the ghost with the
  owner's value, and BCListSPC only lists owned points. So that face stays free unless the
  neighbor cell is fixed too.
- Why it matters: the right/back/top faces of a fixed region that touch a processor
  boundary are not constrained; flow can enter the fixed cells. Results depend on the
  number of CPUs.
- Fix idea: loop over owned faces and check both neighbor cells (needs ghosted phase ratio
  / flag), or add the ghost contributions with ADD/MAX before the scatter.
- Status: fixed (group C, 1f57d0ed): BCApplyPhase and BCApplyCells replaced by BCApplyFixCells.
  The fixed-cell flag is ghosted (one DA_CEN exchange); each rank fixes the faces it owns when a
  cell on either side is fixed. With incompressible flow one lost face is harmless (mass balance
  keeps it zero); a cell losing two faces, at a corner of the fixed region on a rank-boundary
  corner, let flow through (t05 setup, fixed block x < 0.25, z < 0.5, 8 cores: permeability
  4.80e-2 vs 4.64e-2 on 1 core; now equal). t05 unchanged; fix_cell has no test.

### S-045 [major, confirmed] With `exy_num_periods`, outflow markers are deleted, not wrapped
- Tags: major, confirmed, logic, design
- Where: advect.cpp:278, 1086-1096; bc.cpp:906 (commit 6124aa7d)
- What: `actx->periodic` is set for `ExyNumPeriods`, but bc.cpp:906 forbids a periodic grid
  together with Exy. So the x-neighbours are not periodic, `FDSTAGGetPointRanks` gives
  grank = -1, and the marker is deleted after its x-shift.
- Why it matters: in simple-shear runs markers flowing out are lost and not refilled unless
  marker control is on. The flag also forces BASIC_EULER (advect.cpp:280) for no gain.
- Fix idea: really wrap (periodic neighbours in the exchange), or drop `ExyNumPeriods` from the
  flag and document that outflow markers are deleted.
- Status: fixed (group C): periodic marker advection is separate from the periodic mechanical BC
  (Anton). With Exy the marker neighbour list is now periodic in x, so outflow markers wrap to
  the opposite side (also on one CPU). Marker count stays constant (t18, 30 steps, no marker
  control: 13824 kept; before 336 lost). t18 unchanged.

### S-046 [major, possible] Solve and marker advection use different dt when dt is not fixed
- Tags: major, possible, physics, design
- Where: LaMEMLib.cpp:715; tssolve.cpp:308; advect.cpp:661,775,801,811,855,929 (commit 6124aa7d)
- What: the SNES solve (elastic `I2Gdt`, temperature residual) uses the current `ts->dt`.
  `ADVSelectTimeStep` then sets `ts->dt = ts->dt_next` (when `fix_dt` is off) before
  `ADVAdvect`, so marker positions, APS/ATS increments, stress rotation and time use the new dt.
- Why it matters: stress and temperature increments belong to a step of length dt_old but are
  placed at dt_new. dt grows 10% per step by default and can drop a lot on a CFL cut, so
  viscoelastic loading rates are off. Pure viscous flow is not affected. Mechanism confirmed;
  size of the effect not measured.
- Fix idea: advect with the solve dt and apply `dt_next` at the next step (as `fix_dt` does).
- Status: fixed (group C): elastic models already used fixed dt (`fix_dt`), so the claim
  about viscoelastic loading was wrong. The mismatch was in models with temperature but no
  elasticity (temperature solved with the old dt). `fix_dt` now also covers `act_temp_diff`
  (solve depends on dt: elasticity, temperature, kinematic blocks; FSSA not counted, only a
  stabilizer): advect with the solve dt, restart above CFLMAX. Purely viscous models keep the
  exact CFL cut. Results of tests with temperature and no elasticity changed (t12 error vs
  analytic 0.0336 -> 0.0243).
### S-047 [major, confirmed] Surface observations (`Ap = 3`) use wrong grid positions and crash with several CPUs in z
- Tags: major, confirmed, parallel, memory, design
- Where: src/adjoint.cpp:2505-2690; FindPointInCellAdjoint at 102-112
- What:
  - Loops run over corner (DA_COR) indices but write into face arrays (DA_X/DA_Y/DA_Z) with
    z weights from node coordinates. Vx and Vy live at cell centres in z (and Vx at cell centres
    in y, Vy in x), so the weights sit half a cell off. The corner velocity `vgrid` is computed
    and never used.
  - The "point in this CPU" test uses the global z range `gcrdbeg/gcrdend`. A CPU that doesn't
    hold the surface calls FindPointInCellAdjoint with a non-local z. That function uses
    `SETERRQ` but returns `PetscInt`, so the error code (83) comes back as a cell index and is
    used to write into the array.
- Why it matters: wrong results with one CPU in z; out-of-bounds writes with more.
- Fix idea: use local `crdbeg/crdend`, make FindPointInCellAdjoint report errors properly, and
  build weights on each face grid with its own coordinates.
- Status: open

### S-048 [major, possible] Parameter vectors have 100 entries per CPU
- Tags: major, possible, parallel, design
- Where: src/adjoint.cpp:1152-1157 (P, Lb, Ub, grad, val), 1186-1189 (aop->vx..sty)
- What: `VecCreateMPI(PETSC_COMM_WORLD, _MAX_PAR_, ...)` sets the local size, so the global size
  is `100 * nproc`, with the parameters copied on every CPU plus unused zero entries. TAO works on
  this whole vector: norms and dot products grow with CPU count, so `-tao_gatol` etc. mean
  different things on 1 vs. N CPUs.
- Why it matters: inversion stopping and steps may depend on the number of CPUs. Not run in
  parallel to measure it.
- Fix idea: local size `mdN` on rank 0 and 0 elsewhere (or a sequential vector), and broadcast
  the values before setting options.
- Status: open

## Minor

### M-001 [minor, confirmed] Adiabatic heating efficiency is applied twice
- Tags: minor, confirmed, logic, quick
- Where: src/JacResTemp.cpp:525 and :534
- What: `Ha = AdiabHeat*(T*alpha*...)`, then `Ha = AdiabHeat*Ha`. The output value
  (`svBulk->Ha`) has the factor once, but the energy equation gets it squared.
- Why it matters: no effect for 0 or 1. With `Adiabatic_Heat = 0.5`, the model uses 25%.
- Fix idea: remove line 534.
- Status: open

### M-002 [minor, confirmed] Top pressure shift isn't area-weighted and overwrites `p_shift`
- Tags: minor, confirmed, physics, quick
- Where: src/JacRes.cpp:558-602
- What: the top cell layer is averaged with equal weights. With `act_p_shift` on, the result
  replaces any user-given `p_shift` without a warning.
- Why it matters: on grids with variable spacing, the mean top pressure isn't really zero.
  The user's input is silently ignored.
- Status: open

### M-003 [minor, confirmed] `Phase_Melt` density ignores the `mfmax` cap
- Tags: minor, confirmed, logic, quick
- Where: src/constEq.cpp:734-742
- What: in the `Phase_Diagram_melt` branch, `mf` is capped by `mfmax` but then not used.
  The density uses `Pd->mf` directly.
- Fix idea: use `mf` on line 741.
- Status: open

### M-004 [minor, confirmed] Healing defaults set only when printing is on
- Tags: minor, confirmed, logic, quick
- Where: src/phase.cpp:223-232
- What: `APSheal2 = APS2` and `healTau2 = healTau` are set inside `if(PrintOutput)`. The
  adjoint re-reads the database with `PETSC_FALSE` (adjoint.cpp:2881), so that copy has no
  second healing stage.
- Fix idea: move the default assignment out of the print block.
- Status: open

### M-005 [minor, confirmed] Sediment model 3 assumes the box starts at x = 0
- Tags: minor, confirmed, logic, quick
- Where: src/surf.cpp:1531
- What: `dz_x = -dz/BoxWidth*x + dz1` uses x, not x - bx. The rate is rate1 at the left edge
  and rate2 at the right edge only if bx = 0. The comment above it (Gemmer Gaussian) is
  copied from model 2 and is wrong: this is a linear ramp. Unlike model 1, there's no level check.
- Fix idea: use (x - bx).
- Status: fixed (group A)

### M-006 [minor, confirmed] Initial surface perturbation options have unit and NaN problems
- Tags: minor, confirmed, input, parallel, quick
- Where: src/surf.cpp:1711-1754
- What: `-FreeSurf_Wavelength` and `-FreeSurf_AmplNoise` are not scaled (only AmplCos is
  divided by the length scale), so they're wrong unless the length scale is 1. With noise only
  (wavelength 0), 2*pi/0 gives inf, cos(inf) = NaN, and 0*NaN = NaN in every node. Each rank uses
  the same seed, so the noise pattern repeats on every CPU tile and depends on CPU count.
- Fix idea: scale all three, skip the cosine if wavelength is 0, and seed by global node index.
- Status: open

### M-007 [minor, confirmed] Topographic diffusion boundary is not zero-flux
- Tags: minor, confirmed, physics, quick
- Where: src/surf.cpp:1566-1567, 1651, 1662
- What: boundary nodes get no x (or y) term, but their neighbours still exchange with them. That
  isn't a Neumann condition (which would need a half-cell term at the boundary node), and the
  total volume isn't kept. The comment says "zero-flux Neumann".
- Fix idea: use a ghost node with z[-1] = z[1], or fix the comment.
- Status: open

### M-008 [minor, confirmed] Erosion and sedimentation levels overshoot by up to one step
- Tags: minor, confirmed, logic, quick
- Where: src/surf.cpp:1004-1008, 1068-1072, 1346-1350
- What: a node just above the erosion level is lowered by the full rate*dt, so it can end up to
  rate*dt below the level. Sedimentation does the same thing in the other direction.
- Fix idea: clamp to the level, e.g. z = max(z - dz, level).
- Status: fixed (group A)

### M-009 [minor, possible] Surface outside the box silently gets zero velocity
- Tags: minor, possible, logic, quick
- Where: src/surf.cpp:461 (`z >= bz && z < ez`)
- What: a surface node at or above the top of the box (or below the bottom) isn't found on any
  rank, so its velocity is 0. FreeSurfAdvectTopo doesn't clamp, and FastScape clamps to exactly
  the top (fastscape.cpp:1957, 2093), which is excluded by `z < ez`.
- Fix idea: clamp z into [bz, ez] on the top/bottom ranks before the lookup, or warn.
- Status: fixed (group C, low value)

### M-010 [minor, confirmed] Memory leaks in FastScape code
- Tags: minor, confirmed, memory, quick
- Where: src/fastscape.cpp:1210, 1234-1242 vs 1310-1317; 3100-3110
- What: VecScatterCreateToZero makes the output Vecs on every rank, but they're destroyed only on
  rank 0, so each step leaks 3 Vecs per non-root rank (4 at step 0). FastScapeDestroy frees only
  one of `ncoor_ori_x` / `ncoor_ori_y`, but 3D non-uniform grids allocate both.
- Status: open

### M-011 [minor, possible] Non-uniform FastScape interpolation can skip a node
- Tags: minor, possible, logic, quick
- Where: src/fastscape.cpp:1055-1079
- What: if FindIndexForInterpolation misses a point (e.g. the last node, a bit outside the range
  due to roundoff), that node isn't written, so it keeps its old value (0 at the start).
- Fix idea: clamp coordinates into the range, or report an error.
- Status: open

### M-012 [minor, confirmed] `interpStretch` can leave stretch factors unset after the last control polygon
- Tags: minor, confirmed, logic, memory, quick
- Where: src/tools.cpp:628 (used by marker.cpp:1696)
- What: the check is `CtrlPoly[last] != numCtrlPoly-1` but should be `!= numPoly-1`. If the last
  control polygon happens to sit at index numCtrlPoly-1 (e.g. control polygons 1 and 2 of a
  10-slice volume), slices after it get no Sx/Sy. Those are uninitialized stack values that
  stretchPolygon then uses.
- Fix idea: compare with numPoly-1, or simply always fill the tail.
- Status: fixed (group A)

### M-013 [minor, confirmed] Phase 0 can't be used as permeable / plume mantle inflow phase
- Tags: minor, confirmed, input, quick
- Where: src/bc.cpp:848
- What: the check `!bc->phase_inflow_bot` treats phase 0 as "not set", because the default
  is 0 (zeroed struct).
- Why it matters: models where the mantle is phase 0 stop with a wrong error.
- Fix idea: default to -1 and test `< 0`.
- Status: open

### M-014 [minor, confirmed] `relax_dist = 0` with `bvel_face_out` gives NaN at top/bot
- Tags: minor, confirmed, logic, quick
- Where: src/bc.cpp:2097-2099, 2161-2163
- What: with the default `relax_dist = 0`, a cell center with `z == top` (or `z == bot`)
  passes the ramp test and computes `velin - (velin/0)*0 = NaN`.
- Why it matters: rare, but grids with round coordinates can hit it exactly.
- Fix idea: skip the ramp when `relax_dist <= 0`.
- Status: open

### M-015 [minor, confirmed] Plume temperature: different Gaussian width; unset Tbot used as -1
- Tags: minor, confirmed, physics, quick
- Where: src/bc.cpp:1519 vs 3017; 1480-1519; 3005-3017
- What: in 3D, the bottom temperature BC uses `exp(-r²/(2R²))`, while injected markers get
  `exp(-r²/R²)`. 2D uses `R²` in both. Also, the plume profile is built on `Tbot` even when
  no `temp_bot` is given (`Tbot = -1`), as long as `temp_top` is set; injected markers do
  the same.
- Why it matters: injected markers and the BC disagree on temperature in 3D; nonsense
  temperatures if `temp_bot` is missing.
- Fix idea: use one formula; require `temp_bot` when `Plume_InflowBoundary` is on.
- Status: open

### M-016 [minor, confirmed] Typos in BC options are silently ignored; `exy` drops no-slip on vx
- Tags: minor, confirmed, input, quick
- Where: src/bc.cpp:675-679, 712-713, 1784-1802
- What: an unknown `bvel_face` (e.g. "left") leaves `face = 0`, so no inflow BC is applied
  and no message. Same for `bvel_temperature_inflow`. With `exy_strain_rates` set, the
  `else` branch that applies no-slip on vx at front/back/bottom/top is skipped, so those
  no-slip flags are silently ignored for vx.
- Fix idea: error on unknown strings; warn or error when `noslip` and `exy` are combined.
- Status: fixed (group B, batch 1)

### M-017 [minor, confirmed] Coarse edge viscosity at upper boundaries averages in zero ghost cells
- Tags: minor, confirmed, logic, quick
- Where: src/matData.cpp:497-498, 529-530, 561-562 (MatDataRestrictParam3D)
- What: `mnx = tnods - 1` equals the number of fine cells, and fine cell indices run
  0..mnx-1. At the last coarse node `I = 2*i = mnx`, so `if(I > mnx) I--` never fires and the
  code reads ghost cell `mnx`, which is 0 (zeroed in MatDataInitParam, not filled by
  LOCAL_TO_LOCAL). The edge viscosity on top/right/back boundary edges is about half the true
  value (a quarter at corner edges). The lower side clamps correctly. With periodic x the
  lower clamp is also wrong (it should wrap).
- Why it matters: only used by rediscretized coarse levels (`-gmg_mat_free_levels > 0`).
  Gives a worse coarse operator near the top boundary; matters with an open top, where the
  shear term at the top edge does not cancel.
- Fix idea: `if(I >= mnx) I--;` (same for J, K), matching the clamp in matBFBT.cpp:141.
- Status: fixed (group A)

### M-018 [minor, confirmed] Some SNES/KSP failure reasons are not printed or not treated as severe
- Tags: minor, confirmed, logic, quick
- Where: src/nlsolve.cpp:430-515 (SNESPrintConvergedReason)
- What: `SNES_DIVERGED_DTOL`, `SNES_DIVERGED_TR_DELTA` and `SNES_DIVERGED_JACOBIAN_DOMAIN`
  print the "FAILED" header but no reason. Under `SNES_DIVERGED_LINEAR_SOLVE`,
  `KSP_DIVERGED_PC_FAILED` (e.g. zero pivot in LU, failed coarse solve) is not in the severe
  list, so with `-snes_continue_on_fail` the run goes on with a broken preconditioner.
- Fix idea: add the missing cases; treat `KSP_DIVERGED_PC_FAILED` as severe.
- Status: fixed (group A)

### M-019 [minor, possible] MFFD mat-vecs count as SNES function calls
- Tags: minor, possible, logic, quick
- Where: src/nlsolve.cpp:284
- What: the MFFD operator calls `SNESComputeFunction`, which increments the SNES function
  counter. Every Krylov iteration in a Newton step adds one call. PETSc stops with
  `SNES_DIVERGED_FUNCTION_COUNT` when `-snes_max_funcs` is reached (default 10000 in the
  PETSc versions I know of; not checked against the installed 3.22/3.24 source).
- Why it matters: long Newton phases with many linear iterations could stop for this reason.
- Fix idea: call `JacResFormResidual` directly in the MFFD function, or raise
  `snes_max_funcs` at setup. Check the PETSc default first.
- Status: won't fix (group C): mitigated by `-snes_max_funcs unlimited`

### M-020 [minor, confirmed] ADVCheckCorners reads an empty buffer
- Tags: minor, confirmed, memory, quick
- Where: advect.cpp:1485, 1526-1539 (commit 6124aa7d)
- What: if the cell and all its local neighbours have no markers, `n == 0`, the parent loop
  never runs, and `markers[sind]` is read from a zero-size buffer with an unset `sind`.
- Why it matters: crash or garbage marker in emptied in/outflow regions, more likely at rank
  corners.
- Fix idea: skip the corner (or use a default phase) when n == 0.
- Status: fixed (group A)

### M-021 [minor, confirmed] AVD never refills empty volumes
- Tags: minor, confirmed, logic, quick
- Where: AVD.cpp:753-754, 1136-1137 (commit 6124aa7d)
- What: `ninj += min(nmin-n, n)`: with n = 0 nothing is injected, and at most n per pass.
- Why it matters: under `mctrl = avd` an empty volume stays empty, with no warning.
- Status: not an issue (AVD batch): no parent marker for an empty volume; sparse ones fill over steps

### M-022 [minor, confirmed] `FindPointInCell` returns an error code as a cell index
- Tags: minor, confirmed, logic, quick
- Where: AVD.cpp:830-833 (commit 6124aa7d)
- What: `SETERRQ` inside a function returning `PetscInt` returns `PETSC_ERR_USER`, which is
  then used as I/J/K.
- Why it matters: a non-local marker silently corrupts data (`numMarkCell[]`, AVD.cpp:908)
  instead of stopping.
- Fix idea: return the index through an argument; return `PetscErrorCode`.
- Status: fixed (group A)

### M-023 [minor, confirmed] RK2 deletes markers whose midpoint leaves the domain
- Tags: minor, confirmed, logic, quick
- Where: cvi.cpp:204-205 (commit 6124aa7d)
- What: `ADVelDeleteOutflow` runs on the midpoint position; a midpoint exactly at the upper
  wall also counts as outside (`x >= ex`).
- Why it matters: markers near walls with large dt are lost even if the final position is inside.
- Fix idea: clamp the midpoint to the domain.
- Status: fixed (group C, low value)

### M-024 [minor, confirmed] `ADVGetSedPhase` uses stale cell numbers for passive tracers
- Tags: minor, confirmed, logic, quick
- Where: subgrid.cpp:673; LaMEMLib.cpp:731-734; passive_tracer.cpp:843 (commit 6124aa7d)
- What: called after `ADVExchange` but before `ADVMapMarkToCells`, so `actx->cellnum` no longer
  matches the markers.
- Why it matters: the sedimentation phase for passive tracers comes from wrong cells.
- Status: open

### M-025 [minor, confirmed] Temperature-file setup loops over the wrong marker count
- Tags: minor, confirmed, memory, quick
- Where: marker.cpp:604-609, 620-642 (commit 6124aa7d)
- What: `ADVMarkSetTempFile` loops over `nmarkx*nmarky*nmarkz`, not `actx->nummark`. With
  `msetup = files` the counts can differ (write past nummark/markcap, or markers left unset).
  `Ix+1` reads out of range when xp == ex.
- Fix idea: loop to `actx->nummark`; clamp indices to [0, n-2].
- Status: fixed (group A)

### M-026 [minor, confirmed] `ADVMarkSave` crashes when the save path has no '/'
- Tags: minor, confirmed, memory, quick
- Where: marker.cpp:294 (commit 6124aa7d)
- What: `(*strrchr(path, '/')) = '\0'` dereferences NULL for `mark_save_file = mdb`. The file
  is written with header -1, so APS is dropped. (Possible: restart with `msetup = files` may
  fail the exact-markers-per-cell check.)
- Status: partly fixed (group A): APS still not saved

### M-027 [minor, confirmed] Periodic x-shift uses `>` while the rank test uses `>=`
- Tags: minor, confirmed, logic, quick
- Where: advect.cpp:1089 vs fdstag.h:327 (commit 6124aa7d)
- What: a marker with X[0] == ex is sent to the wrap neighbour but not shifted; that rank
  aborts with "Non-local point". Needs exact equality, so rare.
- Fix idea: `X[0] >= ex`.
- Status: fixed (group A)

### M-028 [minor, confirmed] Adiabatic gradient uses |z - Ztop|
- Tags: minor, confirmed, physics, quick
- Where: advect.cpp:2053; bc.cpp:2952 (commit 6124aa7d)
- What: `Z = PetscAbs(P->X[2]-Z_Top)`; rock markers above the reference level (topography)
  are heated as if they were below it.
- Fix idea: `Z = PetscMax(Z_Top - z, 0)`.
- Status: open

### M-029 [minor, confirmed] Ridge spreading rate can silently be 0
- Tags: minor, confirmed, input, quick
- Where: marker.cpp:1193-1200 (commit 6124aa7d)
- What: without `v_spread`, `bc->velin` is used, which is 0 without `bvel_velin` or with
  `velin_array`. The age becomes inf/huge, so the whole box gets T_top.
- Fix idea: require v_spread > 0.
- Status: open

### M-030 [minor, confirmed] Marker input not checked
- Tags: minor, confirmed, input, quick
- Where: parsing.cpp (`getIntParam` only checks max); marker.cpp:378, 536; advect.cpp:156-180;
  cvi.cpp:1239-1240 (commit 6124aa7d)
- What: negative phase IDs pass (out-of-range read at marker.cpp:536, caught later in
  `ADVCheckMarkPhases`); `nmark_sub <= 0`, `nmin > nmax`, `nmark_avd <= 0`, bad
  `nmark_lim` are accepted; `stagp_a = 0` is silently replaced by 2/3.
- Status: fixed (group B, batch 1)

### M-031 [minor, confirmed] Memory leak in `ADVelReAllocStorage`
- Tags: minor, confirmed, memory, quick
- Where: cvi.cpp:800, 826-828 (commit 6124aa7d)
- What: `vi->cellnum` is freed but `vi->markind` is overwritten without `PetscFree`.
- Status: fixed (group A)

### M-032 [minor, confirmed] Dikes without free surface read `phRat[-1]`
- Tags: minor, confirmed, memory, quick
- Where: src/dike.cpp:352 (`GetDikeContr`), 673 (`Compute_sxx_magP`); surf.cpp:37
  (commit 6124aa7d)
- What: `AirPhase` is -1 when `surf_use = 0`, and is used as an index without a check.
- Why it matters: out-of-bounds read; `dikeRHS` and the stress integral get a random extra
  term.
- Fix idea: guard with `AirPhase >= 0`, or require a free surface when `act_dike = 1`.
- Status: open

### M-033 [minor, possible] Dynamic dike can produce NaN
- Tags: minor, possible, logic, quick
- Where: src/dike.cpp:1095, 753-754 (commit 6124aa7d)
- What: the local dike azimuth is `atan(sumslope/sumadd)`. If no other dike cell lies
  within `filty` (cell size in y larger than `filty`, default 1.5 km; or a one-row dike),
  `sumadd = 0` and `azim = NaN`, which spreads to all smoothed stresses. In
  `Compute_sxx_magP`, a column with no rock below `Tsol` has `liththick = 0` → divide by 0.
- Fix idea: fall back to `azim = 0` when `sumadd == 0`; guard `liththick > 0`.
- Status: open

### M-034 [minor, confirmed] Restart with dikes re-reads time stepping
- Tags: minor, confirmed, restart, quick
- Where: src/dike.cpp:1588, 1593, 219 (commit 6124aa7d)
- What: `DynamicDike_ReadRestart` calls `TSSolCreate` again, so if `dt` is in the input the
  current time step is reset to the initial `dt`, and the dt schedule is rebuilt. This
  happens only when `act_dike = 1`. `DBDikeCreate` also resets `istep_count`, so after
  restart new values go into a different slot of the averaging history.
- Fix idea: drop the `TSSolCreate` call (ts is already in the restart file); keep
  `istep_count` from the restored struct.
- Status: open

### M-035 [minor, confirmed] Dike opening rate and units
- Tags: minor, confirmed, input, quick
- Where: src/dike.cpp:306-345 (and 414-453), 212, 855, 1230 (commit 6124aa7d)
- What: dike divergence uses `|bc->velin|` as spreading rate. Without `bvel_velin`
  (e.g. spreading set by other BCs) it is 0 and the dike silently does nothing.
  `magPwidth` is not divided by `scal->length` (`filtx`, `filty` are).
- Fix idea: make the spreading rate an explicit dike input or check it is > 0; scale
  `magPwidth`.
- Status: open

### M-036 [minor, confirmed] Tracer melt phase search checks the wrong marker
- Tags: minor, confirmed, logic, memory, quick
- Where: src/passive_tracer.cpp:599, also 370, 607, 1049 (commit 6124aa7d)
- What: in the loop over markers of the cell, the phase test uses `actx->markers[ii].phase`
  (loop counter) instead of `actx->markers[id_m].phase`. If no marker passes the test, `dist` is
  empty and `dist.begin()->second` reads past the end. The same empty-list read happens in
  `ADV_Assign_Phase` and `Check_advection_condition` if a cell has no markers.
- Fix idea: use `id_m`; check `dist.size()` before use.
- Status: open

### M-037 [minor, confirmed] CFL and CFLMAX range checks can never trigger
- Tags: minor, confirmed, input, quick
- Where: src/tssolve.cpp:56, 61 (commit 6124aa7d)
- What: `if(ts->CFL < 0.0 && ts->CFL > 1.0)` is never true; should be `||`. Same for CFLMAX.
  CFL > 1 or negative values are accepted.
- Status: fixed (group A)

### M-038 [minor, confirmed] Automatic coarse-grid CPU count uses the local cell count
- Tags: minor, confirmed, parallel, quick
- Where: src/options.cpp:482 (commit 6124aa7d); fdstag.cpp:1645 (local sizes)
- What: with `coarse_num_cpu = -1` the target is `ceil(coarse_num_local_cells / coarse_cells_per_cpu)`.
  `coarse_num_local_cells` is the count on one rank, so the result is almost always 1 and the
  whole coarse grid goes to one CPU, however large it is.
- Fix idea: multiply by the number of ranks (total coarse cells).
- Status: fixed (group A)

### M-039 [minor, confirmed] `out_ptr` overwrites the `Passive_Tracer` switch after tracers are set up
- Tags: minor, confirmed, logic, quick
- Where: src/paraViewOutPassiveTracers.cpp:32 (commit 6124aa7d)
- What: PVPtrCreate reads `out_ptr` straight into `jr->ctrl.Passive_Tracer`. It runs after
  `ADVPtrPassive_Tracer_create` (LaMEMLib.cpp:192 vs 204), which uses the same flag.
  - `out_ptr = 1` with `Passive_Tracer = 0` (the default): the flag becomes 1, but no tracer
    vectors were made. Tracer advection (passive_tracer.cpp:465) and output then use NULL Vecs,
    so the run stops with a PETSc error. info/options/input_file.dat sets both keys in different
    places (440, 676), so this is easy to hit.
  - `Passive_Tracer = 1` with `out_ptr = 0`: the flag becomes 0, so tracers are made but silently
    never advected and not written to the restart file.
- Why it matters: an output switch changes the model, or crashes it.
- Fix idea: give PVPtr its own `outptr` flag (the field exists and is unused) and only write
  output when both flags are on.
- Status: open

### M-040 [minor, confirmed] `UpdatePVDFile` uses the file name after `free` in its error message
- Tags: minor, confirmed, memory, quick
- Where: src/paraViewOutBin.cpp:770-772 (commit 6124aa7d)
- What: `free(fname)` comes before `if(fp == NULL) SETERRQ(..., "cannot open file %s", fname)`.
- Why it matters: when the .pvd file can't be opened (e.g. restart in a folder without the old
  .pvd file, which opens it with "r+b"), the error message reads freed memory, which can print
  garbage or crash.
- Fix idea: move `free(fname)` after the check (as in the other open calls).
- Status: fixed (group A)

### M-041 [minor, possible] After a restart from an older checkpoint, .pvd files keep old entries after `</VTKFile>`
- Tags: minor, possible, restart, quick
- Where: src/paraViewOutBin.cpp:766-797 (commit 6124aa7d)
- What: the .pvd offset is stored in the restart file (whole LaMEMLib struct). On restart the
  file is opened with "r+b", the new entry and closing tags are written at the saved offset, but
  the file isn't truncated. If the first run went past the checkpoint, the old later entries
  stay after `</VTKFile>` until the new run writes past them.
- Why it matters: the .pvd file isn't valid XML for a while (or for good, if the new run is
  shorter). ParaView may refuse it or ignore the tail; not tested.
- Fix idea: truncate the file after writing the closing tags (`ftruncate(fileno(fp), ftell(fp))`).
- Status: fixed (group A)

### M-042 [minor, confirmed] `Adjoint_ScaleCostFunction = Var` uses only the last point
- Tags: minor, confirmed, logic, quick
- Where: src/adjoint.cpp:920-947
- What: `var = pow(...)` instead of `+=`. The mean uses `abs(sum)`, so the variance is taken
  around the wrong value for negative data. `mdI = 0` (Ap = 2/3) or zero mean divides by zero.
  "Mean" divides the squared misfit by the mean (not its square) in non-dimensional units, so
  the result still carries units.
- Status: open

### M-043 [minor, confirmed] Scaling-law exponents are wrong for `log10` parameters
- Tags: minor, confirmed, physics, quick
- Where: src/adjoint.cpp:3024-3040, 3080-3095, 3151
- What: `b = grad * P / F` assumes `grad = dF/dP`. For `log10 = 1`, `P` is log10(p) and `grad`
  is `dF/dlog10(p)`, so `b` should be `grad / (F ln 10)`, and `A` should use `p`, not `log10(p)`.
  "Velocity check" is always equal to F by construction. The file writes non-dimensional
  coordinates.
- Status: open

### M-044 [minor, confirmed] A parameter value of 0 gives a zero FD step
- Tags: minor, confirmed, logic, quick
- Where: src/adjoint.cpp:1866 (brute-force FD), 1995 (adjoint dr/dp), 1328 (descent start step)
- What: steps are relative (`eps * value`). A parameter equal to 0 (e.g. a `log10` value of 0,
  or cohesion 0) divides by zero. A zero gradient gives an infinite start step in the built-in
  descent.
- Fix idea: use `eps * max(|value|, 1)` or an absolute step.
- Status: open

### M-045 [minor, confirmed] Dike strain-rate correction is missing in ghost cells
- Tags: minor, confirmed, parallel, medium
- Where: src/JacRes.cpp:1004-1017
- What: the dike part is removed from `dxx/dyy/dzz` in place, only for owned cells, inside
  the cell loop. The ghosts were already exchanged in JacResGetEffStrainRate. The edge
  loops read neighbour cells through ghosts, so at a processor boundary they see
  uncorrected values.
- Why it matters: the edge strain-rate invariant (and so viscosity and plasticity) next to
  processor boundaries differs from the serial run. With dikes, results depend on the
  number of CPUs.
- Fix idea: apply the correction in JacResGetEffStrainRate, before LOCAL_TO_LOCAL. Store
  `dikeRHS` per cell, since gres needs it too.
- Status: open

### M-046 [minor, possible] Dike opening still enters stabilization stress and shear heating
- Tags: minor, possible, physics, medium
- Where: src/constEq.cpp:798-822 (called from JacRes.cpp:1079)
- What: cellConstEq uses `svCell->dxx`, which is not dike-corrected, for the stabilization
  stress `2*eta_st*dxx` and for the dissipative strain rate in `Hr`. The corrected value is
  used only for the main stress.
- Why it matters: dike opening creates false stabilization stress and false shear heating
  in dike cells.
- Status: fixed (group C, 2026-10-27): cellConstEq removes the dike opening (from dikeRHS) from
  the total strain rate used for the stabilization stress and Hr. Changes t26 (only |Div|_inf
  failed in dyndike_4core and dike_variableM; all t26 references updated)

### M-047 [minor, possible] `init_lith_pres` only sets the pressure history; weak stop test
- Tags: minor, possible, logic, medium
- Where: src/JacRes.cpp:1545-1703
- What: the computed lithostatic pressure is copied to `svBulk.pn` and to the markers, but
  not to the pressure in `gsol`, which stays 0. The convergence test compares global sums
  of |p|. Local changes can cancel out. If the sum is 0 (no gravity), you get 0/0 = NaN,
  it runs all iterations and prints a false warning.
- Why it matters: the option name suggests the solution pressure is initialized, but it
  isn't. Check whether that's intended.
- Status: fixed (group C, 92b20ac3): setting only the pressure history is intended (Anton); zero
  pressure (no gravity) now counts as converged, no NaN or false warning

### M-048 [minor, confirmed] Phase temperature `T = 0` (°C) is treated as "not set"
- Tags: minor, confirmed, input, medium
- Where: src/phase.cpp:798-803; src/marker.cpp:523-536
- What: `T`, `T_Nu`, `T_liq` and `T_sol` are shifted to Kelvin only when they are non-zero,
  and 0 means "unset". A user who sets air temperature `T = 0` [°C], which is common, gets
  no phase temperature at all. The code already warns for the same problem with `TRef_fk`
  (phase.cpp:703).
- Fix idea: use a DBL_MAX sentinel for "unset", like `gwLevel`.
- Status: not an issue (Anton, 2026-10-27): air temperature is not set through phase `T`;
  air markers get `Ttop` (advect.cpp, marker.cpp). Phase `T` is meant for hot bodies (plumes,
  intrusions); a phase block at exactly 0 [C] is rare and no input uses it. No change.

### M-049 [minor, possible] Moduli averaged with non-elastic phases counted as 0
- Tags: minor, possible, physics, medium
- Where: src/constEq.cpp:603-626 (G), :697 and :756 (Kb)
- What: `Gavg = sum(phRat*G)` counts phases without G as 0. In a cell that's 50% elastic rock
  and 50% non-elastic phase (often air), I2Gdt doubles. The phase viscosity, though, uses
  each phase's own G. Kb works the same way: a mixed cell becomes more compressible than
  either phase.
- Why it matters: elastic stress history is too strong in mixed cells, for example at the
  free surface. It's hard to see, but it isn't consistent.
- Status: fixed (group C, 92b20ac3): G averaging kept (history stress averaged the same way);
  Kb now averages compressibility 1/Kb (non-elastic = incompressible). Changes t10 (2eb1e72b)

### M-050 [minor, confirmed] Slope erosion and topographic diffusion ignore periodic x
- Tags: minor, confirmed, logic, medium
- Where: src/surf.cpp:1167, 1223, 1643
- What: `I1 < 0 -> I1 = I` also runs for periodic grids, so the left edge is one-sided (no
  x-diffusion there), while the right edge wraps through the ghost node. FreeSurfAdvectTopo
  handles `periodic` correctly; these two functions don't.
- Status: open

### M-051 [minor, confirmed] Rank-0-only errors can hang; restart reads aren't checked
- Tags: minor, confirmed, parallel, restart, medium
- Where: src/tools.cpp:293-296 (DirMake), 319-322, 342-345; tools.cpp:80, fdstag.cpp:292-294
- What: DirMake/DirRemove/DirRename raise SETERRQ on rank 0 only. In DirMake the other ranks
  are waiting in MPI_Barrier, so a failed mkdir hangs the run instead of stopping it. Restart
  readers ignore the return value of fread, so a truncated restart file loads garbage without
  any message.
- Fix idea: broadcast the status and fail on all ranks; check fread counts.
- Status: fixed (batch 4): no hang, PETSc aborts all ranks; fread counts checked

### M-052 [minor, confirmed] Sign of `velin` on Right/Back faces flips depending on `bvel_face_out`
- Tags: minor, confirmed, input, medium
- Where: src/bc.cpp:2103-2123 (x), 2167-2187 (y)
- What: with `bvel_face_out` set, Right/Back use `-vel` (positive `velin` = inflow). Without
  it, they use `+vel` (positive `velin` = outflow). The mass-balance `velout` uses the same
  formula in both cases.
- Why it matters: the same input means opposite flow directions; easy to set up wrong.
- Fix idea: use one convention (inflow positive) in both branches, and document it.
- Status: open

### M-053 [minor, confirmed] Fixed-cell file read: unchecked I/O, one-rank errors hang
- Tags: minor, confirmed, parallel, input, medium
- Where: src/bc.cpp:1229-1274 (BCReadFixCell), 1148 (BCReadRestart)
- What: each rank opens its own file; `SETERRQ(PETSC_COMM_WORLD, ...)` on one rank (missing
  file, wrong size) can hang the others. `stat` and `fread` results are not checked; `fp`
  is not closed on the size error. Same pattern as M-051.
- Fix idea: collect the status with an MPI reduction and fail on all ranks; check counts.
- Status: fixed (batch 4): no hang, PETSc aborts all ranks; stat/fread checked, file closed

### M-054 [minor, possible] CompensatingInflow has no mass-balance check
- Tags: minor, possible, input, medium
- Where: src/bc.cpp:2132-2142, 2213-2223
- What: face 5 sets `vx = velin` on the left and `-velin + 2*velin_net` on the right inside
  [bot, top]. Net inflow is `2*velin*(top-bot)`, which must leave through `bvel_velbot` /
  `bvel_veltop` (default 0) or an open top. Nothing checks this. Also `velin_net` acts as a
  uniform shift of both sides, not as a net inflow as the name suggests.
- Why it matters: a closed box with net inflow has no incompressible solution.
- Fix idea: check the balance at setup, or compute `velbot` from it; document `velin_net`.
- Status: fixed (group C, low value)

### M-055 [minor, confirmed] Matrix-free operator treats the periodic x seam as a free-slip wall
- Tags: minor, confirmed, logic, medium
- Where: src/matFree.cpp:637-640, 683-684 (operator), 1304-1309, 1355-1360 (diagonal);
  src/matBFBT.cpp:208-209 (commit 6124aa7d)
- What: `SET_VEL_TPC` and the diagonal's ghost flags fire at `i == 0` / `i == mnx` with no
  `periodic` check. At the seam this drops the ghost term and sets the wall factor to 0, so
  dvy/dx and dvz/dx are 0 there (free slip). JacResFormResidual and the assembled matrix
  handle periodic correctly. Separately, wBFBT ignores x-velocity constraints when
  `periodic` is on (`cfv[0]`, `cfv[1]` stay 1), so internal fixed-velocity regions leak.
- Why it matters: with `periodic = 1` plus `use_mat_free_jac` (`-js_mat_free`) or
  `-gmg_mat_free_levels`, the Picard operator or MG smoother is wrong at the seam. The
  residual is still right, so a converged result is right, but convergence slows or fails.
  Nothing blocks this combination, and no test covers it.
- Fix idea: skip the x-direction checks when `fs->periodic` is set (as JacRes does), or
  refuse periodic + matrix-free at setup.
- Status: fixed (group B, batch 5)

### M-056 [minor, confirmed] Rediscretized coarse levels lose no-slip and pressure BCs on boundary ghosts
- Tags: minor, confirmed, logic, medium
- Where: src/matData.cpp:642-645 (3D), 808-811 (2D)
- What: coarse BC vectors are set to DBL_MAX and only owned points are filled from the fine
  grid. Ghost points outside the domain (where the fine grid stores no-slip and open-boundary
  pressure values) stay DBL_MAX. On coarse levels that rebuild the operator (assembled level
  below matrix-free levels, and matrix-free coarse levels), `getTwoPointConstr` and
  `SET_VEL_TPC` then see free slip, and `SET_PRES_TPC` sees no pressure BC.
- Why it matters: with `-gmg_mat_free_levels`, the coarse operator differs from the fine one
  at no-slip walls and open boundaries; MG convergence suffers. Galerkin levels are not
  affected (they use R*A*P).
- Fix idea: copy the fine ghost constraint flags to the coarse ghosts (e.g. check the fine
  ghost at the matching boundary position).
- Status: fixed. 3D in group B (batch 5); 2D 2026-10-27, same copy of boundary ghost flags (no effect today: 2D has no rediscretized levels, R/P don't read these ghosts)

### M-057 [minor, confirmed] ADVCheckCorners depends on CPU count; same seed; centre-plane gap
- Tags: minor, confirmed, parallel, medium
- Where: advect.cpp:1368-1376, 1412-1413, 1459, 1568 (commit 6124aa7d)
- What: the parent is searched in local cells only; a new `PetscRandom` with the default seed
  is made every call on every rank; a marker with xp == xc counts in no corner, so markers
  are injected into non-empty corners; the count is printed for rank 0 only.
- Fix idea: use ghost markers or document; seed once per rank; use `>=` on one side; reduce
  the count.
- Status: fixed (group C, low value)

### M-058 [minor, confirmed] AVD half-volume rule uses processor boundaries
- Tags: minor, confirmed, parallel, medium
- Where: AVD.cpp:747-749, 791-793 (commit 6124aa7d)
- What: `(i == 0) | (i+1 == mv->M)` uses the local count, so a face volume cut by a processor
  boundary is checked as two halves needing nmin/2+1 each. `nmax` is not halved for half
  volumes.
- Why it matters: marker control depends on CPU count; boundary volumes can hold twice the
  density before deletion.
- Status: fixed (AVD batch)

### M-059 [minor, confirmed] AVD "edge" control volumes are face volumes; labels swapped
- Tags: minor, confirmed, logic, medium
- Where: AVD.cpp:944-946, 865-884, 813-816 (commit 6124aa7d)
- What: for each dir only one coordinate is shifted to cell centres, so the volume is centred
  on a face, not an edge. The log prints dir 0 as "XYED" and dir 2 as "YZED", the reverse of
  lines 695-697.
- Why it matters: volumes around shear-stress points are never checked; log is misleading.
- Status: fixed (AVD batch)

### M-060 [minor, confirmed] Polygon file sizes not checked against buffers
- Tags: minor, confirmed, memory, medium
- Where: marker.cpp:1466, 1622-1624, 1664-1744 (commit 6124aa7d)
- What: `Polys[_max_polygons_]` (10) overflows with more polygons in one slice; `Vol.num`,
  `len`, `Fcount` are not checked against the allocated sizes; VLAs `SxAll[Vol.num]`.
- Why it matters: a bad or old geomIO file overwrites memory.
- Status: fixed (group B, batch 2)

### M-061 [minor, confirmed] Corner output at the periodic x seam uses one side only
- Tags: minor, confirmed, logic, medium
- Where: interpolate.cpp:102, 158, 220, 377 (commit 6124aa7d)
- What: in periodic mode `i == mx` never happens, so only the `i == 0` clamp fires and the
  corner value at x = 0 ignores ghost cell -1. Output only; velocity (`use_bound = 1`) is fine.
- Status: fixed (x clamps skipped when periodic; ghost -1 holds the wrapped value)

### M-062 [minor, confirmed] Halfspace temperatures use hard-coded kappa = 1e-6
- Tags: minor, confirmed, physics, medium
- Where: marker.cpp:1013, 1152, 1215; phase_transition.cpp:1197, 1339 (Box halfspace T)
  (commit 6124aa7d)
- What: ignores phase conductivity and heat capacity; in non-dimensional mode it becomes a
  non-dimensional 1e-6.
- Fix idea: make kappa an input, or reject in non-dimensional setups.
- Status: open

### M-063 [minor, confirmed] Rank-local hard errors can hang; sedimentation-phase error on fast uplift
- Tags: minor, confirmed, parallel, medium
- Where: subgrid.cpp:85, 606-608; advect.cpp:1906; tools.cpp:293, 798 (commit 6124aa7d)
- What: "Empty control volume", "wrong phase" and `DirMake` errors are raised on one rank in
  collective code, so other ranks hang. subgrid.cpp:606 aborts with "Incorrect sedimentation
  phase" when the surface rises more than one cell per step or has overhangs.
- Status: fixed: hang part not real (PETSc aborts all ranks); sedimentation-phase stop kept (Anton) with a clearer message; `-2` mark no longer overwrites overhang rock cells

### M-064 [minor, confirmed] `ResetParam = APS` is applied every step and inconsistently
- Tags: minor, confirmed, logic, medium
- Where: src/phase_transition.cpp:883-905, 922, 1026-1029 (commit 6124aa7d)
- What:
  - The reset runs every step for every marker that is on the reset side, not once when the
    marker changes phase. So in the new phase APS can never grow.
  - `InsideAbove = 1` means "above" (new phase) for Constant/Clapeyron but "inside" (PH1,
    the "below" phase) for Box. The docs say "reset on particles below PT or within box";
    for Constant it resets above.
  - `Check_NotInAirBox_Phase_Transition` never sets `InAbove`, so it is always 0: with
    `BothWays`/`InsideToOutside` the reset never happens; with `OutsideToInside` APS is
    reset for all listed phases everywhere, every step. The `PhaseOutside = -1` +
    `OutsideToInside` feature (line 922) also never fires for NotInAirBox.
- Why it matters: plastic strain is wiped out where it shouldn't be, or not reset where it
  should.
- Fix idea: reset only when `P->phase` actually changes; return `InAbove` from the
  NotInAirBox check.
- Status: open

### M-065 [minor, confirmed] `BoxVicinity` uses wrong or stale bounds
- Tags: minor, confirmed, logic, medium
- Where: src/phase_transition.cpp:851-864, 393-414 (commit 6124aa7d)
- What: `dzBox` is computed from `bounds[3]-bounds[2]` (y size) instead of
  `bounds[5]-bounds[4]`. For NotInAirBox, `bounds[]` is only set when `nsegs` is not given;
  with `nsegs` it stays 0, so with `BoxVicinity = 1` no marker is ever changed. Moving,
  linked or dynamic-dike boxes also keep their initial `bounds[]`. The docs say "*2 in all
  directions", the code adds one box width on each side.
- Why it matters: wrong set of markers is updated near the box; in thin 2D models the z
  margin can be almost zero. Option is off by default.
- Fix idea: use the z size; for NotInAirBox derive the vicinity from the current
  segment bounds.
- Status: open

### M-066 [minor, confirmed] Time-step schedule can overflow and cuts nstep_max
- Tags: minor, confirmed, memory, medium
- Where: src/tssolve.cpp:381-423 (TSSolMakeSchedule), 457-462 (TSSolAdjustSchedule) (commit 6124aa7d)
- What: `schedule`, `steps` and `ts->schedule` hold `_max_num_steps_` = 2000 values, but the
  number of steps from the periods is never checked. A schedule with more than 2000 steps writes
  past the arrays (`ts->schedule` is inside TSSol, so CFL, time etc. get overwritten).
  `TSSolAdjustSchedule` starts the shift loop at `min(maxSteps, 1999)` and writes
  `schedule[i+1]`, which is index 2000 (out of bounds). `nstep_max` is also set to the schedule
  length, so the run stops at the last period time even if `time_end` is later, without a message.
  The comment says the first `time_dt_periods` entry must be 0, but this is not checked.
- Fix idea: check the total step count against `_max_num_steps_` (or allocate dynamically); keep
  the last `step_dt_periods` value after the schedule ends instead of cutting `nstep_max`.
- Status: fixed (group B, batch 2)

### M-067 [minor, confirmed] Time-step restart repeats phase transitions
- Tags: minor, confirmed, logic, medium
- Where: src/LaMEMLib.cpp:713 (`if(restart) continue;`), src/phase_transition.cpp:967-968 (commit 6124aa7d)
- What: when a fixed time step exceeds CFLMAX, the loop starts the step again from the top.
  `Phase_Transition` runs again, so `MovingBox` moves the box a second time (old dt, then new dt),
  and dike time averaging counts the step twice. With adjoint modes, the gradient for the
  rejected step was already computed and added.
- Fix idea: do the CFL check before applying step-start changes, or skip them on a restart.
- Status: open

### M-068 [minor, confirmed] Passive tracers: APS frozen, wrong IDs, silent fallbacks
- Tags: minor, confirmed, logic, input, medium
- Where: src/passive_tracer.cpp (commit 6124aa7d)
- What:
  - APS is copied from the closest marker at start (376) and never updated (the `aps` array in
    `ADVAdvectPassiveTracer` is unused), but it is written to output (`out_ptr_APS`).
  - ID is `i + ny*j + ny*nx*k` (262); should be `i + nx*j + nx*ny*k`. IDs repeat when nx != ny.
  - An unknown `PassiveTracer_ActiveType` (typo) silently means "Always" (struct is zeroed).
    `Time` exists in the code but not in the docs; `_Phase_ptr_` can't be selected.
    `PassiveTracer_ActiveValue = 0` is rejected, so "T > 0 C" or "time > 0" can't be used.
  - A box partly outside the domain is not checked; those tracers get -DBL_MAX for good. A tracer
    that lands exactly on the upper domain bound (clamp uses `>`, ownership uses `<`) is lost too.
- Status: open

### M-069 [minor, confirmed] SHmax/EHmax direction vectors are averaged to corners, which breaks near north-south
- Tags: minor, confirmed, physics, medium
- Where: src/outFunct.cpp:742-743, 767-768; src/JacResAux.cpp (JacResGetSHmax/EHmax sign rule)
  (commit 6124aa7d)
- What: each cell stores a unit direction with the sign rule `x >= 0`. Directions are axial
  (n and -n mean the same), but they are then averaged linearly from 8 cells to a corner. Near
  north-south, neighbours like (0.01, 1) and (0.01, -1) average to about (0.01, 0), which points
  east-west and is very short.
- Why it matters: wrong stress/extension directions in the output where SHmax is close to the y
  axis. Output only.
- Fix idea: average the doubled-angle form (cos 2θ, sin 2θ) and convert back, or write the field
  as cell data.
- Status: fixed (group B, batch 3)

### M-070 [minor, confirmed] `AllMaterialParameters` scan: writes past arrays, checks only the first block
- Tags: minor, confirmed, memory, medium
- Where: src/adjoint.cpp:204-260, 296-389
- What: the loop over `<AdjointParameterStart>` blocks never does `fb->blockID++`, so only block
  0 is checked for `AllMaterialParameters` (and its `ExcludePhase` lines are read `nblocks`
  times). The count `*iP` is checked against `_MAX_PAR_` only after all writes to
  `type_name/phsar/Par/...` (fixed arrays of 100), so many phases overflow the stack arrays.
  `ExcludePhase` has no `[` check and no limit. `par_str` is `PetscMalloc`ed per profile and
  never freed. All values are parsed with `strtod`, also non-numeric ones.
- Fix idea: advance blockID; check the limit inside AddParamToList.
- Status: open

### M-071 [minor, confirmed] Adjoint solve leaves the forward KSP on the `as_` prefix and cells in a perturbed state
- Tags: minor, confirmed, logic, medium
- Where: src/adjoint.cpp:1935, 1947, 1995-2024
- What: the adjoint uses the SNES's own KSP and sets its prefix to `as_`. It is never set back,
  and the adjoint runs every time step, so all later forward solves read `as_` options. After the
  parameter loop the material DB is restored, but `svCell/svEdge` still hold values from the
  last perturbed residual (`FormResidual`); these go to output and the next step (tiny effect,
  eps = 1e-6).
- Fix idea: use a separate KSP for the adjoint; recompute the residual after the loop.
- Status: open

### M-072 [minor, confirmed] `SyntheticForwardRun` saves a zero vector; the `OFdef = 0` path doesn't exist
- Tags: minor, confirmed, cleanup, medium
- Where: src/adjoint.cpp:1123-1132, 62-66, 2742; LaMEMLib.cpp:690
- What: in `SyntheticForwardRun` the adjoint objective is never called, so `IOparam->xini` stays
  zero and `Forward_Solution_Vel.bin` is all zeros (it would only hold point targets anyway).
  The header describes loading `Forward_Solution.bin` with `OFdef = 0`, but nothing loads it,
  and `Adjoint_ObjectiveFunctionDef` can only be 1.
- Fix idea: save `jr->gsol` and add the load, or remove the mode and the comment.
- Status: open

### M-073 [minor, confirmed] Dynamic dike parallel limits are not checked up front
- Tags: minor, confirmed, parallel, design
- Where: src/dike.cpp:111-120, 959-1001, 1314, 1395 (commit 6124aa7d)
- What: cpu_x = 1 is only checked in `Set_dike_zones`, after `Smooth_sxx_eff` has already
  indexed `DA_CELL_1D` as if it owned all x. Neighbour messages use the local `ny`, so all
  y-ranks must have the same ny. Smoothing only looks one rank up and down in y, so a
  filter wider than a rank gives CPU-dependent results. `DA_CELL_2D_tave` is sized with the
  first dynamic dike's `istep_nave`; another dike with a different value stops with
  "Problems: istep_nave". The new DMDAs use PETSc's default partition, not the FDSTAG one.
- Fix idea: check all limits in `DBDikeCreate`, pass FDSTAG ownership ranges to the DMDAs.
- Status: open

## Notes

### N-001 [note, confirmed] `dikeRHS` is not set when dikes are off
- Tags: note, confirmed, logic, quick
- Where: src/JacRes.cpp:917, passed at :1079
- What: it's only read inside cellConstEq when `actDike` is on, so it's safe today.
  Compilers warn about it, and any change to the `if` order would read garbage.
- Fix idea: initialize it to 0 at the declaration.
- Status: open

### N-002 [note, confirmed] Material data changed inside compute loops; Controls copied per cell
- Tags: note, confirmed, cleanup, quick
- Where: src/JacResAux.cpp:389-390 (`mat->rp` clamped), src/JacResTemp.cpp:107-110
  (`M->nu_k = 1`), src/JacResTemp.cpp:85 (`Controls ctrl = jr->ctrl` copied on every call)
- Why it matters: input checks belong in setup. The struct copy runs for every cell on every
  residual evaluation, which costs time for no reason.
- Status: fixed (group C): rp clamp removed (already checked in phase.cpp); nu_k default no longer written to the material; Controls by pointer

### N-003 [note, confirmed] `eta_min` is not a cutoff
- Tags: note, confirmed, physics, quick
- Where: src/phase.cpp:1511, src/JacRes.cpp:150, src/constEq.cpp:375, 798-800, 894
- What: `eta_min` is only used as the default stabilization viscosity `eta_st`. That's added
  in parallel (`2*eta_st*D`) to every phase's stress and viscosity. It's printed as
  "Minimum viscosity".
- Why it matters: the effective viscosity is eta + eta_min, not max(eta, eta_min). Stress
  can be above the yield stress by 2*eta_min*DII. Fine if eta_min is small; misleading
  otherwise. The docs should say so.
- Status: won't fix (group C): soft eta_min/eta_max limits are intended (Anton)

### N-004 [note, confirmed] Creep input checks never fire; `Bd`/`Bn` silently overwritten
- Tags: note, confirmed, input, quick
- Where: src/phase.cpp:526-548
- What: the diffusion check lists all four eta/Bd combinations, so it can never fail. If both
  `eta` and `Bd` are given, `Bd` is replaced by `1/(2 eta)`. The dislocation check allows
  `eta0+e0+n` together with `Bn`, and `Bn` is then overwritten. The error messages claim
  the parameters must be unique.
- Fix idea: reject (eta && Bd) and (eta0 && Bn).
- Status: fixed (group C, 2c3f19d3): phases are cleared before reading, so a restart takes
  them from the input file only (re-reading was added to override parameters, Anton); checks
  reject eta + Bd and eta0/e0/n + Bn. Default eta_st moved to DBMatSetDefaultEtaSt (JacResCreate
  and restart).

### N-005 [note, confirmed] meltParam.cpp is dead code
- Tags: note, confirmed, cleanup, quick
- Where: src/meltParam.cpp (whole file)
- What: no function in it is called. The only reference is a commented-out line in
  constEq.cpp:144-147.
- Fix idea: remove it or wire it in, but don't leave it half-done.
- Status: open

### N-006 [note, confirmed] PrintMatProp prints wrong values under some labels
- Tags: note, confirmed, logic, quick
- Where: src/phase.cpp:1471, 1479
- What: it prints `rho_c` twice (`rho_n` is missing). In the plasticity line, `eta_st` is
  printed under the "ch" label and `eta_vp`/`eta_st` are swapped. Only used for debugging.
- Status: open

### N-007 [note, confirmed] Small clean-ups in surf/fastscape
- Tags: note, confirmed, cleanup, quick
- Where / what:
  - surf.cpp:111-112 reads `marginE` twice; surf.cpp:244-249 the same block twice.
  - surf.cpp:72-73 and 1275: comment says "non-dimensional value of 1 m/yr", but `unit_m_yr`
    is the m/yr value of one non-dimensional unit (code is right).
  - fastscape.cpp:92 checks `<= 2` but the message says "must be >= 2".
  - fastscape.cpp:129-130 and fastscape.h:163: kd is in m^2/yr, not m/yr.
  - fastscape.cpp:124 says kf, kd "can be set as an array", but only one value is read.
  - fastscape.cpp:1528 builds the output folder time as t*s + dt*s, while LaMEM uses (t+dt)*s
    (LaMEMLib.cpp:542). A rounding difference can give two folders for one step.
- Status: fixed (group C): all items; FastScape folder time now (t+dt)*s

### N-008 [note, confirmed] Average topography is a plain node average
- Tags: note, confirmed, physics, quick
- Where: src/surf.cpp:288-292
- What: `avg_topo` is the mean over nodes, not weighted by area. On non-uniform grids,
  infinitely fast erosion (model 1) doesn't keep the volume.
- Status: fixed (group C, low value)

### N-009 [note, confirmed] ISA / eigenvalue code is unused and has latent bugs
- Tags: note, confirmed, cleanup, quick
- Where: src/Tensor.cpp:112-318 (Tensor2RNEigen, getISA), 493-603 (Tensor2RSSpectral)
- What: none of these functions has a caller (`_ISA_` exists only as a name in objFunct.cpp).
  Tensor2RNEigen treats |D| < tol as "three zero eigenvalues", but D = 0 also for repeated
  non-zero eigenvalues (uniaxial flow, e.g. diag(2,-1,-1)), which would be reported as simple
  shear. Tensor2RSSpectral gives NaN for an all-zero tensor (max < 0 is never true, so it
  divides 0/0).
- Fix idea: remove, or fix before using it again (test p and q, not D; return early when nrm = 0).
- Status: fixed (group C): removed

### N-010 [note, confirmed] Non-unit bias is rejected, so all bias code is dead
- Tags: note, confirmed, cleanup, quick
- Where: src/fdstag.cpp:80-83
- What: MeshSeg1DReadParam stops with "Non-unit bias ratios are deprecated". So the biased
  branches in MeshSeg1DGenCoord (fdstag.cpp:153-171), FSLoadNonUniformGrid and
  FastScapeCreateGlobalGrid (fastscape.cpp:578-592, 921-937) never run. `bias_x` is still
  read as an input.
- Fix idea: remove the dead branches, or drop the check if bias should come back.
- Status: fixed (group C): bias removed fully; `bias_*` keys no longer read

### N-011 [note, confirmed] Velocity cylinder/box edge cases
- Tags: note, confirmed, input, quick
- Where: src/bc.cpp:2405-2411, 2234, 2369
- What: with `vmag` and base == cap the axis length is 0, giving NaN velocities; radius 0
  divides by zero. No checks. Velocity boxes and cylinders are skipped during the initial
  guess solve, but Bezier blocks and boundary inflow are not.
- Status: fixed (group C): radius > 0 and base != cap checked; initial-guess skip of boxes/cylinders left as is

### N-012 [note, confirmed] Wrong comments and messages, debug print, small clean-ups
- Tags: note, confirmed, cleanup, quick
- Where: bc.h:263 and others
- What:
  - bc.h:263 says `Gaussian=0=default, Poiseuille=1`; code uses Poiseuille=0, Gaussian=1.
    `Plume_VelocityType` is `_REQUIRED_` but passes a default "Gaussian" (unused).
  - Error messages: "Influx_type" should be "Inflow_Type" (782); dimension error says
    "Plume_Type" (827).
  - BCGetVelins prints a debug line every step when both period arrays are used, with a
    hard-coded "(Myr)" label (1683).
  - "Bezier" blocks are linear in time between path points; `path` array is sized for
    3n-2 Bezier points.
  - `// nondimensionalize` comment on `velin_net` (2070) is wrong; `velin_net` lookup is
    duplicated from BCGetVelins.
  - `k_thermal = 1e-6` m²/s is hard-coded in BCOverridePhase (2961).
  - Summary prints plume velocity type for the permeable plume, where it is not set; the
    periodic message is printed in the middle of the no-slip mask line.
  - `face_out = -1` (inflow on both sides): the opposite face gets no phase/temperature
    override in BCOverridePhase.
  - `tSPCList/tSPCVals` are allocated but never used; `vel = 0.0` duplicated (2158).
  - TPC shear uses `SIZE_CELL(j, sy, ...)` with the ghost-range start (1788). Index is off
    by one on the last rank, but the ghost spacing equals the last cell, so no effect.
- Status: open

### N-013 [note, possible] AVD point with no cells gives NaN coordinates
- Tags: note, possible, logic, quick
- Where: AVD.cpp:1123-1125 (commit 6124aa7d)
- What: division by `hclaim` with no guard; `hclaim = 0` only if a point owns no AVD cells
  (unlikely).
- Status: fixed (group A)

### N-014 [note, confirmed] Docs say Clapeyron parameters overwrite a predefined law
- Tags: note, confirmed, input, quick
- Where: info/options/input_file.dat:743-749; src/phase_transition.cpp:562-580
  (commit 6124aa7d)
- What: the example sets `Name_Clapeyron = Mantle_Transition_660km` together with
  `clapeyron_slope`, `P0_clapeyron`, `T0_clapeyron` marked "overwrite". The code reads
  these only for `user-defined`, so they are ignored.
- Fix idea: read them as optional overrides after `SetClapeyron_Eq`, or fix the docs.
- Status: fixed (group C, low value)

### N-015 [note, confirmed] Dead code: whole ObjFunct module and the pointwise FD kernels
- Tags: note, confirmed, cleanup, quick
- Where: src/objFunct.cpp (all), src/adjoint.cpp:3684-4150
- What: `ObjFunctCreate/CompErr/Destroy` are never called. If revived they have bugs: the
  `MPI_Bcast` of `otUse` uses `MPIU_2INT` (twice the size, overwrites the next fields);
  `infile` is never set but is opened and freed; `-objf_use_x 0` still turns a field on; topo
  error is scaled with the velocity scale; BOUG/ISA/SHMAX errors are summed but never computed.
  In adjoint.cpp, cellConstEqFD / edgeConstEqFD / devConstEqFD / setUpPhaseFD (~470 lines) are
  never called; `AdjointGradientResetParameter` is declared but not defined; `aop->gradfield`,
  `aop->Perturb`, `CurScalst`, `Adjoint_Vecs.val`, `DII_ref`, `mdI<_MAX_OBS_` tests are unused
  or always true.
- Fix idea: delete.
- Status: fixed (group C): objFunct.cpp, FD kernels and unused fields removed

### N-016 [note, confirmed] `gw_level_type = surf` gives a flat water table
- Tags: note, confirmed, physics, medium
- Where: src/JacResAux.cpp:357
- What: it uses `surf->avg_topo`, the mean topography, not the local surface height.
- Why it matters: with real topography, pore pressure is too high under valleys and too low
  under hills. Users may expect the water table to follow the surface.
- Status: not an issue (group C, Anton): flat level at average topography is intended; comment added

### N-017 [note, possible] Creep profile experiment types not checked
- Tags: note, possible, physics, medium
- Where: src/phase.cpp:888-1352
- What: each profile sets `_UniAxial_` or `_SimpleShear_`, and that changes the prefactor
  by up to 3^((n+1)/2)/2 vs 2^(n-1). For example, the HK2003 laws are marked simple shear,
  but most of their data are triaxial. Needs a check against each paper.
- Status: not an issue (group C, Anton)

### N-018 [note, confirmed] Polygon helpers use stack arrays and can divide by zero
- Tags: note, confirmed, memory, medium
- Where: src/tools.cpp:608-609 (interpStretch), 658-662 (findCenterMass); marker.cpp:1676-1677
- What: arrays sized by the number of polygon nodes/slices live on the stack (variable-length
  arrays, not standard C++). Large polygons can overflow the stack. findCenterMass divides by
  the area with no check, so a degenerate polygon gives inf/NaN.
- Status: fixed (group C): stack arrays removed; zero-area polygon uses node mean (marker.cpp arrays were already heap)

### N-019 [note, confirmed] Newton/Picard switch logic
- Tags: note, confirmed, logic, medium
- Where: src/nlsolve.cpp:216-238
- What: switching back to Picard tests `nrm > refRes*rtolNwt`, with `refRes` the residual
  at the first iteration, not at the switch. So a Newton step can grow the residual by up to
  `rtolNwt/rtolPic` (120x with defaults) before falling back. After `maxItNwt` Newton steps
  the code does one Picard step and then switches straight back to Newton (the residual is
  still below `rtolPic*refRes`). Printed label says "MMFD" instead of "MFFD".
- Fix idea: confirm this is intended; otherwise store the residual at the switch and use it
  for the divergence test.
- Status: fixed (group C, 92b20ac3): -snes_newton_maxit is the Newton budget per solve, then
  Picard only; MFFD label fixed; divergence test vs first residual of the step kept

### N-020 [note, confirmed] `jp_pgamma` means different things per solver; inv_eta Schur drops a term
- Tags: note, confirmed, physics, medium
- Where: src/lsolve.cpp:48-90; src/matrix.cpp:316, 1080, 1131, 1147
- What: with `-jp_type bf`, `pgamma > 1` augments the velocity block (`getVelSchur`). With
  `-jp_type user`, it only scales the pressure diagonal to `-1/(pgamma*eta)`; no velocity
  augmentation. For the bf inverse-viscosity Schur preconditioner with penalty, `iS = 1/diag
  = -pgamma*eta`; the usual augmented-Lagrangian estimate is `-(1 + pgamma)*eta`. Small
  difference for large `pgamma`. Also: assembled Picard rows for constrained DOFs are
  identity, matrix-free rows are zero (works only because those residuals are always zero).
- Fix idea: document `jp_pgamma` per solver type; consider `-(1+pgamma)*eta`.
- Status: not an issue (group C, Anton)

### N-021 [note, confirmed] Option overrides, stale comments, dead code, small clean-ups
- Tags: note, confirmed, cleanup, medium
- Where: various
- What:
  - nlsolve.cpp:86-93: `SNESSetFromOptions` runs before the `js_` prefix is set, so
    un-prefixed `-ksp_*` options also reach the Stokes KSP; `PCSetType(pc, PCMAT)` after
    `KSPSetFromOptions` silently overrides `-js_pc_type`; `MATMFFD_WP` (286) overrides
    `-fd_mat_mffd_type`. `SNESGetLineSearch` is called twice (79, 84).
  - multigrid.cpp:283-284: `PCSetFromOptions` before `PCMGSetGalerkin(NONE)`, so
    `-gmg_pc_mg_galerkin` is ignored without a message.
  - multigrid.h:88 declares `MGLevelCreate` with 3 arguments; the definition has 4 (C++
    overload that is never defined).
  - multigrid.cpp:1603-1802: `compareVecs`, `genRandVec`, `VecSetBC`, `TestInterp`,
    `TestInterpBC` are never called.
  - matData.cpp:157-161: comment says ghost points store the negative primary index; code
    stores -1 and uses `pdofidx` instead.
  - lsolve.cpp:82-90: bf-only checks (`bf_vs_type`, `bf_schur_type` with `jp_pgamma`) also
    fire for `jp_type user`.
  - matAux.cpp:124/127 `VecZeroEntries` twice; matAux.cpp:65 second `MatSetFromOptions`.
  - Restriction/prolongation weights assume uniform spacing, and R is not a scaled P^T.
- Status: partly fixed (group C): dead multigrid test functions removed, MGLevelCreate declaration fixed; option order, MFFD type, matData comment, lsolve checks, matAux duplicates fixed 2026-10-27 (1a67766e); R/P weights and coarse viscosity averaging moved to improvements.md (not a bug)

### N-022 [note, confirmed] Dead code, wrong comments and log lines, small clean-ups
- Tags: note, confirmed, cleanup, medium
- Where: various (commit 6124aa7d)
- What:
  - Dead: `AVDLoadPoints`, `AVDInjectDeletePoints`, `AVDExecuteMarkerInjection`
    (AVD.cpp:403-665); `_DISP_` (advect.h:71); test block subgrid.cpp:740-786.
  - Wrong comments: cvi.cpp:13,157 ("RK 4th order"); marker.cpp:973 ("sin", code is -cos).
  - Logs with `PETSC_COMM_WORLD` print rank 0's local counts only (AVD.cpp:818,
    subgrid.cpp:174, advect.cpp inject count).
  - `nbyte` cast to `PetscMPIInt` (advect.cpp:1135, cvi.cpp:681/694) overflows at ~15-20M
    markers per neighbour.
  - Ellipsoid parser lacks `fb->ID = jj` (marker.cpp:1073); unknown `Temperature` strings
    ignored; degenerate shapes not checked; polygon noise amplitude differs from InitCoord.
  - `A->buffer = 1` (AVD.cpp:98) grows arrays one at a time (O(n²)).
  - `FreeSurfGetAirPhaseRatio` called twice per step (advect.cpp:1638, LaMEMLib.cpp:776).
  - Nearest-rock search in `ADVMarkCrossFreeSurf` depends on storage order.
- Status: fixed (group C): dead code removed (3c04c606); ellipsoid `fb->ID` already set; marker
  exchange stops with an error above the MPI count (d3a18147); comments fixed (cvi.cpp, marker.cpp
  cosine); AVD and subgrid logs print totals over all ranks (one MPI_Allreduce each); unknown
  `Temperature` strings stop the run; AVD chains double their capacity (old copy read past the
  array for buffer > 1); polygon noise uses the marker spacing as in ADVMarkInitCoord;
  nearest-rock search uses the phases from before the update. FreeSurfGetAirPhaseRatio kept
  twice (only call with advect = none, second call is a no-op). Also fixed: layer `rand_ampl`
  was scaled by maxPhaseID instead of the length unit. All tests pass, no reference changes.

### N-023 [note, confirmed] Phase transitions and dikes: dead code, names, small clean-ups
- Tags: note, confirmed, cleanup, medium
- Where: src/phase_transition.cpp, src/dike.cpp (commit 6124aa7d)
- What:
  - `Locate_Dike_Zones`: `istep+1 == 0` is never true (dike.cpp:508); the "debugging"
    condition (522-523) is left in, so `Compute_sxx_magP` and `Smooth_sxx_eff` run every
    step; `magma_presence = 1.0 //testing` (748).
  - `x_Mc`, `z_Mc` declared and mentioned in comments but never read (dike.h:45-46).
  - `PhaseInside[0] = -1` is set after `Set_Box_Phase_Transition` prints "Don't set phase",
    so that message never shows (phase_transition.cpp:130 vs 333, 451).
  - `Parameter_transition = Depth` is the z coordinate: `PhaseAbove` applies for
    `z >= value` (shallower). Pressure value is always in Pa (even GEO), time comment
    says [s] but value is in time units.
  - `MovingBox` moves only x bounds, ignores partial overlap with `[t0_box, t1_box]`; a
    linked box with a larger ID lags one step (`LinkNotInAirBoxes`).
  - `DynamicPhTr_WriteRestart2` declared but not defined (phase_transition.h:49).
  - `Set_dike_zones`: `ixmax-1` can be `sx-2` when the dike is next to the x boundary.
  - `fread`/`fwrite` results not checked (phase_transition.cpp:1478-1515).
  - Time averaging divides by `istep_nave` from the start, so the first steps are biased
    toward 0 (history starts at 0).
- Status: open

### N-024 [note, confirmed] Input parsing: silent bad values and small risks
- Tags: note, confirmed, input, medium
- Where: src/parsing.cpp (commit 6124aa7d)
- What:
  - `strtol(ptr, NULL, 0)` (509) reads `010` as 8 and `08` as 0; `1.5` or `1e3` for an integer
    gives 1. `strtod` and `strtol` results are never checked, so `abc` becomes 0 silently.
  - Extra values after the expected `num` are silently ignored (507, 568).
  - Strings stop at the first space, so paths with spaces don't work; an `=` inside a value
    (e.g. a file name) fails the "equal signs must be surrounded by spaces" check.
  - `FBLoad` errors (no `-ParamFile`, file can't be opened) happen on rank 0 only, while other
    ranks wait in `MPI_Bcast` (hang). `fread`/`ftell` results are not checked.
  - `getIntParam`/`getScalarParam`/`getStringParam` use `fb->nblocks` before the `&& fb` check.
  - Block ranges are checked only for begin < end of the same index; nested or overlapping blocks
    are not detected.
  - `PetscOptionsReadFromFile` keeps only the first token after the key, so options with several
    values (separated by spaces) are cut silently.
- Status: fixed (group B): approved parts; parser relaxed

### N-025 [note, confirmed] Output and time-step quirks, labels, docs
- Tags: note, confirmed, cleanup, medium
- Where: src/tssolve.cpp, options.cpp, scaling.cpp, LaMEMLib.cpp, doc/src/man/GettingStarted.md
  (commit 6124aa7d)
- What:
  - `nstep_out <= 0` turns off all output, also `dt_out` and the first step (tssolve.cpp:216).
    `time_out` is reset by any output, so with `nstep_out = 1` (default) `dt_out` does nothing.
  - With `fix_dt`, step 0 uses `dt` and the schedule is applied one step late
    (`dt_next = schedule[istep]` is used in the next step).
  - The last step can overshoot `time_end` and output times are not hit exactly (no dt trimming).
  - The steady-state temperature solver (`its_` prefix) uses the `-ts_ksp_atol_auto` flag and
    updates `ts_ksp_ref_norm` of the transient solver (LaMEMLib.cpp:1010-1042).
  - Options set with `%g` keep only 6 digits (options.cpp:907).
  - `lbl_heat_production` says `[W/m3]` but the unit is power/mass, W/kg (scaling.cpp:188, 246).
    "WRNING! ... scaling is employed" has a typo and no newline (scaling.cpp:123).
  - GettingStarted.md:49-60 says `nstep_out` sets the number of time steps; it is `nstep_max`.
- Status: fixed. Group B: labels/typos, fix_dt schedule. Group C (2026-10-27): nstep_out = 0 means
  output only by dt_out (and first/last step), nstep_out < 0 no output; dt_out output times on a
  fixed grid (multiples of dt_out, not reset by other outputs); dike.cpp no longer divides by
  nstep_out = 0. No dt trimming at time_end / output times (Anton: running past the end is fine)

### N-026 [note, confirmed] Output quirks: names, stale residual, dead vectors, unchecked inputs
- Tags: note, confirmed, cleanup, medium
- Where: src/outFunct.cpp, paraViewOutBin.cpp, paraViewOutAVD.cpp, paraViewOutPassiveTracers.cpp
  (commit 6124aa7d)
- What:
  - `plast_dissip` writes `Hr`, the full shear heating (viscous + plastic, constEq.cpp:822, 909),
    not only plastic dissipation.
  - `energ_res` writes `jr->ge` from before the last temperature solve (nlsolve.cpp:349-353), so
    it is not the remaining residual.
  - `vol_rate`, `vorticity`, `ang_vel_mag` are empty stubs and can't be switched on. If they were
    on, the header would promise data that's never written, so the file would be broken.
  - With groundwater on, `pressure` writes total pressure, the same as `total_pressure`.
    `total_pressure` uses one-sided values at the domain edge (`use_bound = 0`) while `pressure`
    uses the boundary ghost values, so the two differ at the edges even without pore pressure.
    Passive tracer pressure (paraViewOutPassiveTracers.cpp:319) has no `p_shift`, grid pressure has.
  - Phase aggregate `phaseID` isn't checked for negative values, so `phase_mask[-1]` is written
    (outFunct.cpp:131). `sprintf(name, "%s %s")` can overflow `name[260]` with a long aggregate name.
  - `out_avd_ref <= 0` isn't rejected (division by zero, paraViewOutAVD.cpp:365). A rank with no
    markers reads `points[-1]` (paraViewOutAVD.cpp:909).
  - Marker, tracer and AVD files are written one value per `fwrite` (slow for many markers).
- Status: fixed. Group B: energ_res, pShift, AVD, outputs. Group C (2026-10-27): total_pressure
  uses the boundary ghost values like pressure (differed at edges only with a pressure BC, e.g.
  pres_top). fwrite per value: won't fix (stdio buffers the writes; ~1-2 s per 10M markers)

### N-027 [note, confirmed] Adjoint input and small quirks
- Tags: note, confirmed, input, medium
- Where: src/adjoint.cpp (several)
- What:
  - `Adjoint_GradientCalculation` and `Adjoint_ScaleCostFunction` are read as optional, but an
    empty value hits the error branch, so they're required; the defaults are never used (449-458).
  - `Vel_comp` is printed uninitialised for `Ap > 1` (903). `Coordinate` keeps the previous
    block's value if missing (804). `Av` is uninitialised for `Ap = 2/3` with fewer than 3
    observation blocks (2485-2508).
  - The final check `!IOparam->Ap` can never be true (Ap is 1..3) (951).
  - `Avel_num` (size `_MAX_OBS_`) is reduced with `_MAX_PAR_` (same value today) (1768, 1776).
  - AdjointPointInPro runs twice per objective call (1609, 1701).
  - Adjoint_ApplyBCs sets constrained entries of dF to the BC values, while the comment says
    they should be zero; the pressure part reads the lists and does nothing (4155-4187). Harmless
    only if the matrix and dr/dp are zero in those rows/columns. Not checked.
  - CreateModifiedMaterialDatabase stores a pointer to a stack `Scaling` in
    `dbm_modified.scal` (2879), and rebuilds the whole DB twice per parameter per time step.
  - Points on a CPU border may be printed twice (2184, 2199).
- Status: open

### N-028 [note, confirmed] FSSA only stabilizes the diagonal terms
- Tags: note, confirmed, physics, design
- Where: src/JacRes.cpp:1086-1103
- What: each momentum equation gets only its own velocity times its own density gradient
  (for example vz*drho/dz*gz). The full term is (v.grad rho)*g, which would include
  vx*drho/dx*gz.
- Why it matters: lateral density gradients (sloped surfaces, tilted interfaces) are not
  stabilized. The code comment says "lumped", but the user docs don't mention the limit.
- Status: won't fix (group C): Boris: leave

### N-029 [note, confirmed] Rock/air phase change keeps old stress, APS and T
- Tags: note, confirmed, physics, design
- Where: subgrid.cpp:528-606 (commit 6124aa7d)
- What: `ADVMarkCrossFreeSurf` changes only the phase; history is mixed into cut cells for one
  step. There is already a WARNING comment.
- Status: won't fix (group C): Boris: leave
