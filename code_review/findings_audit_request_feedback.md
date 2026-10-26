# Findings audit — request for feedback

This is `findings_audit.md` with two additions per finding: **who wrote the affected lines** (from `git blame`) and **who is asked for feedback**. Boris's remarks (`- Boris: ...`) are unchanged. The original text starts after the "Feedback request" section below.

## Feedback request

Asked to review: **Wenrong Cao, Anton Popov, Garrett Ito, Andrea Piccolo, Yun Luo** and Boris Kaus (for his own code and Georg's). For each finding assigned to you, please reply with one of: *confirmed bug* / *not a bug (explain intent)* / *already fixed (commit)* — and, if you can, a fix or a PR. Each finding below lists the commits that wrote the lines, so you can see what you are being asked about. Being named does not mean the bug is your fault; several were introduced by later refactors, or by people who have since left.

Rule used to pick reviewers: any of the six who authored the affected code; code by people outside the group is routed as you asked (Georg Reuber → Boris; Arne Spang, Adina Püsök, Tobias Baumann → Anton; Jana Schierjott → Garrett). If none of them did, the finding is assigned by code area (FastScape/surf → Wenrong, Yun Luo; dike → Garrett; phase transitions, passive tracers, phase/marker setup → Andrea; everything else → Anton) and is marked as such.

| Reviewer | Findings (by authorship) | Findings (by area only) |
|---|---|---|
| Wenrong Cao | 5: S-032, S-033, M-007, M-050, N-007 | 0: – |
| Anton Popov | 80: S-001, S-002, S-003, S-004, S-006, S-007, S-009, S-013, S-014, S-018, S-019, S-024, S-029, S-035, S-036, S-039, S-044, S-045, S-046, M-005, M-008, M-009, M-012, M-016, M-017, M-018, M-019, M-020, M-021, M-022, M-023, M-025, M-026, M-027, M-030, M-031, M-037, M-038, M-040, M-041, M-046, M-047, M-048, M-049, M-051, M-053, M-055, M-056, M-057, M-058, M-059, M-060, M-061, M-063, M-066, M-067, M-069, N-002, N-003, N-004, N-007, N-008, N-009, N-010, N-011, N-013, N-014, N-015, N-016, N-017, N-018, N-019, N-020, N-021, N-022, N-024, N-025, N-026, N-028, N-029 | 3: S-010, S-011, M-054 |
| Garrett Ito | 21: S-015, S-022, S-023, S-037, M-004, M-029, M-032, M-033, M-034, M-035, M-045, M-046, M-048, M-062, M-064, M-065, M-067, M-073, N-001, N-002, N-023 | 0: – |
| Andrea Piccolo | 22: S-008, S-020, S-021, S-030, S-036, S-038, M-001, M-003, M-013, M-014, M-015, M-016, M-024, M-028, M-036, M-039, M-052, M-068, N-005, N-011, N-012, N-023 | 0: – |
| Yun Luo | 8: S-005, S-034, S-043, M-010, M-011, M-036, N-007, N-010 | 0: – |
| Boris Kaus | 43: S-005, S-008, S-012, S-016, S-017, S-020, S-021, S-025, S-026, S-027, S-028, S-031, S-037, S-039, S-040, S-041, S-042, S-047, S-048, M-002, M-006, M-016, M-039, M-042, M-043, M-044, M-062, M-064, M-065, M-070, M-071, M-072, N-003, N-006, N-012, N-015, N-017, N-019, N-021, N-022, N-025, N-027, N-028 | 0: – |

### For Wenrong Cao

| ID | Where (HEAD) | Title | Why you | Boris |
|---|---|---|---|---|
| S-032 | src/surf.cpp:1229-1253 (centered gradient + erosion); slope scan 1173-1188 | Slope-dependent erosion never erodes peaks and ridges | Implemented FreeSurfAppSlopeErosion with centered-difference slope. |  |
| S-033 | src/surf.cpp:1199-1211 (E_max, dt_cfl = 0.5*d_min/E_max, nsolve) | Slope-erosion sub-step limit isn't a stability limit | Wrote the sub-step rule limiting vertical drop to half the grid spacing. |  |
| M-007 | src/surf.cpp:1566-1567,1649-1669 (FreeSurfAppTopoDiffusion) | Topographic diffusion boundary is not zero-flux | boundary nodes skip x/y flux terms; comment claims zero-flux Neumann |  |
| M-050 | src/surf.cpp:1167,1223,1643 | Slope erosion and topographic diffusion ignore periodic x | I1<0 -> I1=I clamp without periodic check (slope erosion, topo diffusion) |  |
| N-007 | src/surf.cpp:111-112,243-248,72-73,1275; src/fastscape.cpp:92,124-130,1528; src/fastscape.h:163 | Small clean-ups in surf/fastscape | unit_m_yr comment wrong in slope-erosion code |  |

### For Anton Popov

| ID | Where (HEAD) | Title | Why you | Boris |
|---|---|---|---|---|
| S-001 | src/JacRes.cpp:1524 | `init_pres` without `fix_phase` reads `phRat[-1]` | wrote phRat[fixPhase] test without checking fixPhase == -1 | please fix |
| S-002 | src/constEq.cpp:588-592 | Softening length scale `Lm` looks inverted | introduced Lm-based APS1/APS2 *= Le/Lm scaling | Anton, please check (hard to believe that this is wring but could be) |
| S-003 | src/phase.cpp:1037 | Profile "Diabase-Caristan_1982" has Vn = 1 m^3/mol | rewrote/moved profiles into phase.c carrying Vn = 1 (blame target) | should be zero (thats how it is in GeoParams - see https://github.com/JuliaGeodynamics/GeoParams.jl/blob/32dfaedb509ee1258551bf59c12a159e250fdfa1/src/CreepLaw/Data/DislocationCreep.jl#L190) |
| S-004 | src/surf.cpp:1444-1462 | Sediment model 2 uses values left over from the previous node | routed: written by Tobias Baumann | yes could be wrong, please fix |
| S-006 | src/fdstag.cpp:550-551 | Grid stretching moves the global box about 0, not about `bg_ref_point` | added gcrdbeg/gcrdend *= (1+eps) in Discret1DStretch (ref=0 assumed) | Please fix & add a test |
| S-007 | src/fdstag.cpp:726-737 | CFL time step ignores negative velocities on non-uniform grids | wrote idt = v/h for non-uniform grid (no abs value) | fix |
| S-009 | src/bc.cpp:1679-1725 (BCGetVelins), 2047 (BCApplyBoundVel) | User `bvel_velout` is always overwritten by the mass-balance value | original mass-balance velout formula (previously only computed if velout unset) | I would anyways like to have a much more general way to set boundary & internal velocities. Adding a flag seems a good idea; default behavior should still be to compute this in a mass-conservative way |
| S-010 | src/bc.cpp:690-699, 1690-1700 | `velin_num_periods = 1` sets the inflow velocity to zero and ignores net periods | by area (written by Iskander Ibragimov) | please fix |
| S-011 | src/bc.cpp:1704-1718 | After restart, net-only inflow adds the net offset twice | by area (written by Panupong Kongpet) | please fix if you think this makes sense |
| S-013 | src/advect.cpp:1217-1228; src/AVD.cpp:591,1205 | Garbage collection deletes the wrong markers when the AVD delete list is unsorted | garbage-collect loop that assumes ascending idel | please fix even when his is no longer default in thye code |
| S-014 | src/cvi.cpp:1305/1314, 1355/1364, 1405/1414 | STAGP interpolation treats processor boundaries as domain walls | added STAGP interpolation with local IN==0/nx wall tests | seems urgent to fix |
| S-018 | src/tools.cpp:628-635; src/marker.cpp:1685-1697 | Polygon stretch factors after the last control polygon are left unset | routed: written by Arne Spang |  |
| S-019 | src/advect.cpp:225-228; src/marker.cpp:210-262 | `rand_noiseGP` moves markers after cell mapping without remapping | routed: written by Tobias Baumann | Please fix |
| S-024 | src/passive_tracer.cpp:560-561 | Passive tracer pressure and temperature are read half a cell too high | renamed variables in same lines; bug preserved |  |
| S-029 | src/JacResAux.cpp:244-321 (lith), 339-411 (pore); edge averages src/JacRes.cpp:1195-1198, 1302-1305, 1410-1413 | Lithostatic and pore pressure are zero in boundary ghost cells | Wrote lp_lith/lp_pore computation (owned cells only, zero ghosts) and pore-pressure edge averaging. |  |
| S-035 | src/AVD.cpp:1186-1214 (AVDDeletePointsMV; area[i]=tclaimed at 1195) | AVD deletion ignores the Voronoi volumes | Split delete path into AVDDeletePointsMV, copying area=tclaimed without the claim-counting loop. |  |
| S-036 | src/phase_transition.cpp:592-594 (scaling), 667-711 (predefined laws in deg C / MPa/K); scaling.cpp | Clapeyron laws are wrong in SI units | routed: written by Arne Spang, Andrea Piccolo |  |
| S-039 | src/options.cpp:937-960 (setSolverOptions clear -> defaults -> reinsert); defaults building scaling/grid at 120-165, 251-256 | Command-line options are ignored while solver defaults are set | Defaults now create Scaling/FDSTAG and read act_temp_diff from input file, ignoring command line. |  |
| S-044 | src/bc.cpp:2594-2606 (BCApplyPhase), 2647-2659 (BCApplyCells) | `fix_phase` / `fix_cell` lose faces on processor boundaries | Wrote the fixed-phase/fixed-cell loops setting i+1/j+1/k+1 faces. |  |
| S-045 | src/advect.cpp:277-280 (periodic flag with ExyNumPeriods), 1084-1090; src/bc.cpp:904-909 | With `exy_num_periods`, outflow markers are deleted, not wrapped | Set actx->periodic for Exy and forbade periodic grid together with Exy. |  |
| S-046 | src/tssolve.cpp:307-308 (if(!fix_dt) ts->dt = dt_next); call at src/LaMEMLib.cpp:697 | Solve and marker advection use different dt when dt is not fixed | Introduced immediate dt update to dt_next before marker advection when dt is not fixed. |  |
| M-005 | src/surf.cpp:1531 (SedimentModel==3) | Sediment model 3 assumes the box starts at x = 0 | routed: written by Tobias Baumann |  |
| M-008 | src/surf.cpp:1004-1008, 1068-1072, 1346-1350 | Erosion and sedimentation levels overshoot by up to one step | routed: written by Xinxin Wang, Wenrong Cao, Tobias Baumann |  |
| M-009 | src/surf.cpp:461 (FreeSurfAdvectTopo velocity lookup) | Surface outside the box silently gets zero velocity | introduced half-open test z >= bz && z < ez |  |
| M-012 | src/tools.cpp:628 (interpStretch) | `interpStretch` can leave stretch factors unset after the last control polygon | routed: written by Arne Spang |  |
| M-016 | src/bc.cpp:675-679, 712-713, 1775-1802 | Typos in BC options are silently ignored; `exy` drops no-slip on vx | added if(Exy) ... else so no-slip vx BCs skipped when exy set |  |
| M-017 | src/matData.cpp:497-498, 529-530, 561-562 (MatDataRestrictParam3D) | Coarse edge viscosity at upper boundaries averages in zero ghost cells | wrote 'if(I > mnx) I--' edge clamp (should be >=) in rediscretized coarsening |  |
| M-018 | src/nlsolve.cpp:430-515 (SNESPrintConvergedReason) | Some SNES/KSP failure reasons are not printed or not treated as severe | severe-KSP list lacks KSP_DIVERGED_PC_FAILED |  |
| M-019 | src/nlsolve.cpp:284 (MatMFFDSetFunction) ; mitigation src/options.cpp:59 | MFFD mat-vecs count as SNES function calls | MFFD operator calls SNESComputeFunction (counted by SNES) |  |
| M-020 | src/advect.cpp:1539 (markers[sind], ADVCheckCorners) | ADVCheckCorners reads an empty buffer | routed: written by Adina Püsök |  |
| M-021 | src/AVD.cpp:753-754 | AVD never refills empty volumes | wrote 'if((nmin-n)>n) ninj += n' logic injecting at most n markers |  |
| M-022 | src/AVD.cpp:823-850 (FindPointInCell, SETERRQ at 830-833) | `FindPointInCell` returns an error code as a cell index | added SETERRQ 'Non-local marker' inside PetscInt-returning function |  |
| M-023 | src/cvi.cpp:204-205 (ADVelRungeKuttaStep) and ADVelDeleteOutflow 477-520 | RK2 deletes markers whose midpoint leaves the domain | routed: written by Adina Püsök |  |
| M-025 | src/marker.cpp:604-607, 614-642 (ADVMarkSetTempFile) | Temperature-file setup loops over the wrong marker count | reworked nmarkx/y/z count lines and interpolation lines |  |
| M-026 | src/marker.cpp:294 | `ADVMarkSave` crashes when the save path has no '/' | added (*strrchr(path,'/'))='\0' without NULL check |  |
| M-027 | src/advect.cpp:1089 vs src/fdstag.h:327 | Periodic x-shift uses `>` while the rank test uses `>=` | periodic wrap uses X[0] > ex |  |
| M-030 | src/advect.cpp:96-109,176-180; src/cvi.cpp:1238-1240; src/marker.cpp:536 | Marker input not checked | nmark_lim/nmark_avd/nmark_sub read with no range checks |  |
| M-031 | src/cvi.cpp:800,825-828 (ADVelReAllocStorage) | Memory leak in `ADVelReAllocStorage` | routed: written by Adina Püsök |  |
| M-037 | src/tssolve.cpp:56,61 | CFL and CFLMAX range checks can never trigger | wrote CFL/CFLMAX range checks with && instead of // |  |
| M-038 | src/options.cpp:482 (+ fdstag.cpp:1666-1671 local size) | Automatic coarse-grid CPU count uses the local cell count | changed coarse CPU estimate from global cell count to coarse_num_local_cells |  |
| M-040 | src/paraViewOutBin.cpp:770-772 | `UpdatePVDFile` uses the file name after `free` in its error message | moved free(fname) before the fp==NULL error using fname |  |
| M-041 | src/paraViewOutBin.cpp:766-797 | After a restart from an older checkpoint, .pvd files keep old entries after `</VTKFile>` | restart-mode r+b open and offset seek without truncating the .pvd |  |
| M-046 | src/constEq.cpp:799-801,817-824; JacRes.cpp:1011-1016 | Dike opening still enters stabilization stress and shear heating | stabilization stress and Hr use svCell->dxx |  |
| M-047 | src/JacRes.cpp:1545-1700 | `init_lith_pres` only sets the pressure history; weak stop test | norm-sum stop test; copies p only to svBulk.pn and markers |  |
| M-048 | src/phase.cpp:786-791; src/marker.cpp:523-536 | Phase temperature `T = 0` (°C) is treated as "not set" | T!=0 treated as set; Kelvin shift only if nonzero |  |
| M-049 | src/constEq.cpp:620-623 (G), 697,756 (Kb) | Moduli averaged with non-elastic phases counted as 0 | Gavg/Kavg summed with non-elastic phases counted as zero |  |
| M-051 | src/tools.cpp:295,321,344,80; src/fdstag.cpp:292-294 | Rank-0-only errors can hang; restart reads aren't checked | rank-0-only SETERRQ in DirMake/DirRemove/DirRename |  |
| M-053 | src/bc.cpp:1148, 1257-1269 | Fixed-cell file read: unchecked I/O, one-rank errors hang | per-rank fixed cell file read; stat/fread unchecked |  |
| M-054 | src/bc.cpp:2129-2142, 2213-2223 | CompensatingInflow has no mass-balance check | by area (written by Thomas Morrow) |  |
| M-055 | src/matFree.cpp:639-640,683-684,1307-1308,1358-1359; src/matBFBT.cpp:208-209 | Matrix-free operator treats the periodic x seam as a free-slip wall | matrix-free two-point constraints and diagonal flags at i==0/mnx ignore periodic |  |
| M-056 | src/matData.cpp:642-645, 808-811 | Rediscretized coarse levels lose no-slip and pressure BCs on boundary ghosts | coarse BC vectors set to DBL_MAX, ghost constraints not restricted |  |
| M-057 | src/advect.cpp:1368-1376,1412-1413,1459,1568 | ADVCheckCorners depends on CPU count; same seed; centre-plane gap | switched random context to PETSC_COMM_SELF |  |
| M-058 | src/AVD.cpp:747-749,791-793 | AVD half-volume rule uses processor boundaries | half-volume rule uses local mv->M/N/P index bounds |  |
| M-059 | src/AVD.cpp:944-946,865-870,814-816 vs 695-697 | AVD "edge" control volumes are face volumes; labels swapped | face-centred volumes called edges; XYED/YZED labels swapped |  |
| M-060 | src/marker.cpp:1466,1643,1676-1677,1730-1732 | Polygon file sizes not checked against buffers | routed: written by Arne Spang |  |
| M-061 | src/interpolate.cpp:102-103,158-159,220-222,377-378 | Corner output at the periodic x seam uses one side only | corner interpolation clamps only at i==0 / i==mx, no periodic handling |  |
| M-063 | src/subgrid.cpp:85,606-608; src/advect.cpp:1906-1908; src/tools.cpp:295,798-800 | Rank-local hard errors can hang; sedimentation-phase error on fast uplift | subgrid errors raised on one rank; sedimentation phase error |  |
| M-066 | src/tssolve.cpp:407-419,457-462 | Time-step schedule can overflow and cuts nstep_max | routed: written by Arne Spang |  |
| M-067 | src/LaMEMLib.cpp:697-700; src/phase_transition.cpp:958-970 | Time-step restart repeats phase transitions | CFL restart "continue" after step-start transitions |  |
| M-069 | src/JacResAux.cpp:77,156; src/outFunct.cpp:742-743,767-768 | SHmax/EHmax direction vectors are averaged to corners, which breaks near north-south | x>=0 sign rule on axial directions; averaged linearly to corners |  |
| N-002 | src/JacResTemp.cpp:67,85,107-110; src/JacResAux.cpp:389-390 | Material data changed inside compute loops; Controls copied per cell | mat->rp clamped inside pore-pressure loop |  |
| N-003 | src/JacRes.cpp:92,150; src/phase.cpp:1499; src/constEq.cpp:375,795-801,891-894 | `eta_min` is not a cutoff | added eta_st stabilization in parallel to all phases; eta_min used as default eta_st (JacRes.cpp:150) |  |
| N-004 | src/phase.cpp:514-523 (diffusion), 527-536 (dislocation) | Creep input checks never fire; `Bd`/`Bn` silently overwritten | check collapsed to always-true (eta&&Bd allowed) in 'formatting' commit; Bd/Bn silently overwritten from eta/eta0,e0 |  |
| N-007 | src/surf.cpp:111-112,243-248,72-73,1275; src/fastscape.cpp:92,124-130,1528; src/fastscape.h:163 | Small clean-ups in surf/fastscape | routed: written by Yun Luo, Wenrong Cao, Tobias Baumann |  |
| N-008 | src/surf.cpp:288-292 (FreeSurfGetAvgTopo) | Average topography is a plain node average | VecSum(gtopo)/(Nx*Ny*Nz) plain node average introduced |  |
| N-009 | src/Tensor.cpp:107-310 (Tensor2RNEigen, getISA), 488-603 (Tensor2RSSpectral) | ISA / eigenvalue code is unused and has latent bugs | ISA/eigen algorithm after Kaminsky 2004; ISA debugged |  |
| N-010 | src/fdstag.cpp:80-83; dead branches fdstag.cpp:144-171, src/fastscape.cpp:564-631,907 | Non-unit bias is rejected, so all bias code is dead | deprecated non-unit bias (error), leaving biased branches dead; original bias mesh code |  |
| N-011 | src/bc.cpp:2399-2411 (BCApplyVelCylinder vmag/a), 2234, 2369 (initGuess skip) | Velocity cylinder/box edge cases | routed: written by Arne Spang, Andrea Piccolo |  |
| N-013 | src/AVD.cpp:1123-1125 (also 546-548) | AVD point with no cells gives NaN coordinates | added (PetscScalar) cast on the division |  |
| N-014 | info/options/input_file.dat:743-746; src/phase_transition.cpp:560-580 | Docs say Clapeyron parameters overwrite a predefined law | made clapeyron_slope/P0/T0 read only for user-defined; overrides dropped |  |
| N-015 | src/objFunct.cpp (all); src/adjoint.cpp:3684-4150 | Dead code: whole ObjFunct module and the pointwise FD kernels | routed: written by Tobias Baumann, Georg Reuber |  |
| N-016 | src/JacResAux.cpp:357 | `gw_level_type = surf` gives a flat water table | gwLevel = surf->avg_topo for gw_level_type surf |  |
| N-017 | src/phase.cpp:840-1352 (SetDiffProfile/SetDislProfile) | Creep profile experiment types not checked | assigned _SimpleShear_/_UniAxial_ types to creep profiles; original creep profile database |  |
| N-018 | src/tools.cpp:608-609,658-662,704-705; src/marker.cpp:1676-1677 | Polygon helpers use stack arrays and can divide by zero | routed: written by Arne Spang |  |
| N-019 | src/nlsolve.cpp:209-232 (switch logic), 247 (MMFD label) | Newton/Picard switch logic | original Newton/Picard switch and 'MMFD' label; maxItNwt/rtolNwt defaults |  |
| N-020 | src/lsolve.cpp:48-90; src/matrix.cpp:316,1080,1131,1147 | `jp_pgamma` means different things per solver; inv_eta Schur drops a term | original pgamma / inv_eta Schur code; jp_type user/bf handling; pgamma checks in lsolve |  |
| N-021 | src/nlsolve.cpp:79-93,286; src/multigrid.cpp:283-284,1603-1802; src/multigrid.h:88; src/matData.cpp:159; src/lsolve.cpp:82-90; src/matAux.cpp:65,124-127 | Option overrides, stale comments, dead code, small clean-ups | SNES/KSP option ordering, PCMAT; unused multigrid test helpers; MGLevelCreate 3-arg declaration |  |
| N-022 | src/AVD.cpp:98,403-665; src/advect.h:71; src/subgrid.cpp:740-786; src/cvi.cpp:13,157; src/marker.cpp:974,1071-1082; advect.cpp:1137; cvi.cpp:683,696 | Dead code, wrong comments and log lines, small clean-ups | subgrid test block; nbyte cast to PetscMPIInt |  |
| N-024 | src/parsing.cpp:370,431 (strtol/strtod), 507-570, FBLoad ~56-87 | Input parsing: silent bad values and small risks | strtol/strtod unchecked parsing; new input framework (FBLoad, rank-0 errors) | many of these make sense and would make the code more robust |
| N-025 | src/tssolve.cpp:216,234-296; src/options.cpp:906-907; src/scaling.cpp:118,183,241; src/LaMEMLib.cpp:986-1026; doc/src/man/GettingStarted.md:49-60 | Output and time-step quirks, labels, docs | its_ steady-state temp solver shares ts_ksp_ref_norm; %g option storing 6 digits; scaling labels / typos | the obvious ones can be fixed; nstep<=0 should indeed switch off output |
| N-026 | src/outFunct.cpp:131,429-488,617-640,672; src/paraViewOutBin.cpp:286; paraViewOutAVD.cpp:365,909; paraViewOutPassiveTracers.cpp:319 | Output quirks: names, stale residual, dead vectors, unchecked inputs | phase aggregate output without negative-ID check; stub output functions; use_bound handling for total_pressure | plast_dissip is correct; energ_res to be fixed; vol_rate, vorticity, ang_vel_mag to get an option to write; Passive tracer pressure should get p_shift; out_avd_ref <= 0 should be rejecetd |
| N-028 | src/JacRes.cpp:1086-1103 | FSSA only stabilizes the diagonal terms | FSSA density-gradient stabilization, diagonal terms only | leave as-is, as this is the key contribution. I'm pretty sure I had an option to do the full terms @ some stage |
| N-029 | src/subgrid.cpp:537-615 (ADVMarkCrossFreeSurf) | Rock/air phase change keeps old stress, APS and T | marker-surface interaction overhaul; phase-only change, WARNING comments; initial marker/free-surface rock-air logic | don't think its easy to figure out what that old stress should be if we go from air->rock (other way is easier), so I would leave as-is |

### For Garrett Ito

| ID | Where (HEAD) | Title | Why you | Boris |
|---|---|---|---|---|
| S-015 | src/marker.cpp:2358 | Oblique ridge axis is wrong unless the front edge is at y = 0 | routed: written by Jana Schierjott | yes, seems wronng; fix |
| S-022 | src/dike.cpp:1331 | Dike stress time averaging loops over the wrong range | x loop bound uses sx+ny instead of sx+nx |  |
| S-023 | src/dike.cpp:463-479 | Dike heating scales the whole cell conductivity by the dike fraction | introduced k = kfac*k with kfac summing only dike phRat |  |
| S-037 | src/phase_transition.cpp:1302-1304 (zbounds[0],[1] test), 1322-1323, 1334; read at 413 | NotInAirBox uses only the first segment's `zbounds` | Made xbounds per segment and zbounds read 2*nsegs but kept using zbounds[0]/[1]. |  |
| M-004 | src/phase.cpp:215-220 (inside if(PrintOutput)) | Healing defaults set only when printing is on | set APSheal2/healTau2 defaults only inside the PrintOutput branch |  |
| M-029 | src/marker.cpp:1193-1200 | Ridge spreading rate can silently be 0 | routed: written by Jana Schierjott |  |
| M-032 | src/dike.cpp:352 (GetDikeContr), 673 (Compute_sxx_magP); src/surf.cpp:37 | Dikes without free surface read `phRat[-1]` | phRat[i]+phRat[AirPhase] with AirPhase possibly -1 |  |
| M-033 | src/dike.cpp:1095 (azim), 753-754 (liththick division) | Dynamic dike can produce NaN | azim = atan(sumslope/sumadd) without sumadd==0 guard; divides by liththick |  |
| M-034 | src/dike.cpp:1588 (TSSolCreate in DynamicDike_ReadRestart), 219 (istep_count reset) | Restart with dikes re-reads time stepping | added TSSolCreate(ts, fb) in restart reader, resetting dt |  |
| M-035 | src/dike.cpp:306,317,340,414,425,448 (v_spread=|bc->velin|); 212,1230 (magPwidth) | Dike opening rate and units | added magPwidth read with scale 1.0 (not length) |  |
| M-045 | src/JacRes.cpp:1008-1016 | Dike strain-rate correction is missing in ghost cells | routed: written by Jana Schierjott |  |
| M-046 | src/constEq.cpp:799-801,817-824; JacRes.cpp:1011-1016 | Dike opening still enters stabilization stress and shear heating | routed: written by Jana Schierjott, Anton Popov |  |
| M-048 | src/phase.cpp:786-791; src/marker.cpp:523-536 | Phase temperature `T = 0` (°C) is treated as "not set" | same non-zero shift for T_liq/T_sol |  |
| M-062 | src/marker.cpp:1014,1152,1215; src/phase_transition.cpp:1197,1339 | Halfspace temperatures use hard-coded kappa = 1e-6 | routed: written by Boris Kaus, Jana Schierjott |  |
| M-064 | src/phase_transition.cpp:886-905,922,1024-1029 | `ResetParam = APS` is applied every step and inconsistently | routed: written by Boris Kaus, Jana Schierjott |  |
| M-065 | src/phase_transition.cpp:851-864 (BoxVicinity), 398-412 (bounds) | `BoxVicinity` uses wrong or stale bounds | bounds[] only set when nsegs not given |  |
| M-067 | src/LaMEMLib.cpp:697-700; src/phase_transition.cpp:958-970 | Time-step restart repeats phase transitions | MovingBox moves celly_xboundL/R by v_box*dt |  |
| M-073 | src/dike.cpp:108-120,959-1001,1314,1397 | Dynamic dike parallel limits are not checked up front | DA_CELL_1D/2D_tave created with default partition, sized from first dike |  |
| N-001 | src/JacRes.cpp:917 (decl), 1008-1011, 1079 | `dikeRHS` is not set when dikes are off | routed: written by Jana Schierjott |  |
| N-002 | src/JacResTemp.cpp:67,85,107-110; src/JacResAux.cpp:389-390 | Material data changed inside compute loops; Controls copied per cell | routed: written by Jana Schierjott, Anton Popov |  |
| N-023 | src/dike.cpp:508,522-523,748,1437; src/dike.h:44-46; src/phase_transition.cpp:130,333,451,975,1478-1515; phase_transition.h:49 | Phase transitions and dikes: dead code, names, small clean-ups | istep+1==0 dead test; debugging condition left in Locate_Dike_Zones; magma_presence=1.0 //testing | Garrett should double check. Pressure is always in Pa, also with GEO units, so that is fine. Most others seem good to be improved by GI should check. |

### For Andrea Piccolo

| ID | Where (HEAD) | Title | Why you | Boris |
|---|---|---|---|---|
| S-008 | src/bc.cpp:723 (read), 1087 (print), 1123 (nondimensionalize) | `bvel_temperature_constant` is scaled with the time unit, not as a temperature | first appearance of bvel_constant_temperature read with scal->time scale | hard to believe that this is wrong.. |
| S-020 | src/phase_transition.cpp:756-763, 811-812, 824-825 | Box with `PhaseOutside = -1` sets marker phase to -1 | introduced nphc flag | Seems indeed a bug; Please fix |
| S-021 | src/phase_transition.cpp:85, 100-119, 227-259, 426-433; dike.cpp:163 | Typos in phase transition input are silently accepted | ID range check and string-compare chains without else error |  |
| S-030 | src/JacResTemp.cpp:517-526 (bdpdx.. lines 520-522 use P[...] incl. ghost neighbours) | Adiabatic heating sees a fake pressure jump at boundaries | Added adiabatic-heating pressure-gradient term using neighbour-cell lithostatic pressure. |  |
| S-036 | src/phase_transition.cpp:592-594 (scaling), 667-711 (predefined laws in deg C / MPa/K); scaling.cpp | Clapeyron laws are wrong in SI units | Added T0 + Tshift conversion and predefined laws with T0 in deg C. |  |
| S-038 | src/LaMEMLib.cpp:720-721 (call before ADVRemap); src/passive_tracer.cpp:585-607, 1032-1049, 904-938 | Passive tracers use a stale marker-to-cell list | Placed ADVAdvectPassiveTracer before ADVRemap and used markstart/markind in tracer code. |  |
| M-001 | src/JacResTemp.cpp:525,534 | Adiabatic heating efficiency is applied twice | added second Ha = AdiabHeat*Ha (line 534) after svBulk->Ha stored |  |
| M-003 | src/constEq.cpp:739-741 (Phase_Diagram_melt branch) | `Phase_Melt` density ignores the `mfmax` cap | added mf cap copy but kept Pd->mf in melt-density line 741 |  |
| M-013 | src/bc.cpp:848 | Phase 0 can't be used as permeable / plume mantle inflow phase | check '!bc->phase_inflow_bot' treats phase 0 as unset |  |
| M-014 | src/bc.cpp:2097-2098, 2161-2162 | `relax_dist = 0` with `bvel_face_out` gives NaN at top/bot | wrote velin-(velin/relax_dist)*(z-top) ramp without relax_dist>0 guard |  |
| M-015 | src/bc.cpp:1507,1519 (BC) vs 3005,3017 (marker injection) | Plume temperature: different Gaussian width; unset Tbot used as -1 | marker T_inflow = Tbot+(Tplume-Tbot)*exp(-r2/R2) with Tbot possibly -1 |  |
| M-016 | src/bc.cpp:675-679, 712-713, 1775-1802 | Typos in BC options are silently ignored; `exy` drops no-slip on vx | bvel_temperature_inflow strcmp chain; typo silently ignored |  |
| M-024 | src/LaMEMLib.cpp:718-721; src/passive_tracer.cpp:843; src/subgrid.cpp:631-ff | `ADVGetSedPhase` uses stale cell numbers for passive tracers | calls ADVGetSedPhase in passive-tracer code using stale actx->cellnum |  |
| M-028 | src/advect.cpp:2053; src/bc.cpp:2952 | Adiabatic gradient uses /z - Ztop/ | introduced PetscAbs(z - Z_Top) in both marker adiabat routines |  |
| M-036 | src/passive_tracer.cpp:599 (ii vs id_m), 370, 606-607, 936, 1048-1049 | Tracer melt phase search checks the wrong marker | mat[actx->markers[ii].phase] uses loop counter; unguarded dist.begin() |  |
| M-039 | src/paraViewOutPassiveTracers.cpp:32 | `out_ptr` overwrites the `Passive_Tracer` switch after tracers are set up | original Passive_Tracer-based PVPtr output switch that out_ptr overwrote |  |
| M-052 | src/bc.cpp:2105-2114 (x), 2169-2179 (y) vs 2122-2123, 2186-2187 | Sign of `velin` on Right/Back faces flips depending on `bvel_face_out` | non-face_out branch (+vel on Right/Back faces) from earlier BC code |  |
| M-068 | src/passive_tracer.cpp:262,376,60-63,523 (aps unused) | Passive tracers: APS frozen, wrong IDs, silent fallbacks | ID formula with ny/nx swapped; silent default Always; ActiveValue==0 rejected |  |
| N-005 | src/meltParam.cpp, src/meltParam.h (whole); src/constEq.cpp:144-147 (commented call) | meltParam.cpp is dead code | commented-out Compute_Melt_Fraction call in constEq.cpp |  |
| N-011 | src/bc.cpp:2399-2411 (BCApplyVelCylinder vmag/a), 2234, 2369 (initGuess skip) | Velocity cylinder/box edge cases | initGuess skip for boxes/cylinders only |  |
| N-012 | src/bc.h:263; src/bc.cpp:782,827,1683,2961 (tSPCList partly used in JacResTemp.cpp:415) | Wrong comments and messages, debug print, small clean-ups | Influx_type message; plume code; hard-coded k_thermal in BCOverridePhase |  |
| N-023 | src/dike.cpp:508,522-523,748,1437; src/dike.h:44-46; src/phase_transition.cpp:130,333,451,975,1478-1515; phase_transition.h:49 | Phase transitions and dikes: dead code, names, small clean-ups | PhaseInside[0]=-1 set after 'Don't set phase' message | Garrett should double check. Pressure is always in Pa, also with GEO units, so that is fine. Most others seem good to be improved by GI should check. |

### For Yun Luo

| ID | Where (HEAD) | Title | Why you | Boris |
|---|---|---|---|---|
| S-005 | src/fastscape.cpp:2776-2804 and 2905-2914 | FastScape restart loses basement and silt fraction if their output is off | copy of basement/silt fraction conditional on output flag; restart always uses them | dont think this was ever tested wiuth restarting; would be good to have |
| S-034 | src/LaMEMLib.cpp:749-752 (FastScapeRun then FreeSurfSmoothMaxAngle); src/fastscape.cpp:1590-1606 (gtopo_fs reused) | With FastScape, max-angle smoothing is undone every step | Added FastScape coupling: smoothing after FastScapeRun while FastScape restarts from own saved topography. |  |
| S-043 | src/fastscape.cpp:456-469 (ranges set once), 646-668 (CreateSurfaceGrid mode 0 never called; only modes 1 and 2 at 178, 2948) | FastScape grid doesn't follow a stretching LaMEM box | Wrote the FastScape grid setup with fixed ranges and unused update-range mode. |  |
| M-010 | src/fastscape.cpp:1210,1234-1242,1305-1317; 3100-3110 (FastScapeDestroy) | Memory leaks in FastScape code | added FastScape coupling file; Vecs from VecScatterCreateToZero leaked off-root; one ncoor_ori freed |  |
| M-011 | src/fastscape.cpp:1046-1081 (if(indices.found) without else) | Non-uniform FastScape interpolation can skip a node | interpolation loop skips node when FindIndexForInterpolation misses |  |
| M-036 | src/passive_tracer.cpp:599 (ii vs id_m), 370, 606-607, 936, 1048-1049 | Tracer melt phase search checks the wrong marker | touched line 599 (pdn test style) but kept markers[ii] |  |
| N-007 | src/surf.cpp:111-112,243-248,72-73,1275; src/fastscape.cpp:92,124-130,1528; src/fastscape.h:163 | Small clean-ups in surf/fastscape | FastScape module: <=2 check, kd units/array comments, time_fs; dup block surf.cpp:243 |  |
| N-010 | src/fdstag.cpp:80-83; dead branches fdstag.cpp:144-171, src/fastscape.cpp:564-631,907 | Non-unit bias is rejected, so all bias code is dead | FastScape grid copies biased branches, unreachable |  |

### For Boris Kaus

| ID | Where (HEAD) | Title | Why you | Boris |
|---|---|---|---|---|
| S-005 | src/fastscape.cpp:2776-2804 and 2905-2914 | FastScape restart loses basement and silt fraction if their output is off | merge resolution turned output copy into out_flags[] array loop (same behaviour) | dont think this was ever tested wiuth restarting; would be good to have |
| S-008 | src/bc.cpp:723 (read), 1087 (print), 1123 (nondimensionalize) | `bvel_temperature_constant` is scaled with the time unit, not as a temperature | added print multiplying by scal->time, masking the wrong scale in the log | hard to believe that this is wrong.. |
| S-012 | src/bc.cpp:3125-3134 | 3D Gaussian plume inflow is not mass-conserving (missing radius² factor) | derived 3D erf integral without the c^2 factor | please fix; may explain some of the issue we had wuth this, |
| S-016 | src/marker.cpp:985, 976, 2207 | Layer `rand_ampl` is divided by the max phase ID | added rand_ampl read with maxPhaseID instead of chLen; rand() noise | not siure how this happened; fix |
| S-017 | src/marker.cpp:1996-2000, 2088-2094, 2043 (and 2054-2086) | Phase diagram load fails silently; size not checked against buffer | wrote T/P range print (minor) | Please fix |
| S-020 | src/phase_transition.cpp:756-763, 811-812, 824-825 | Box with `PhaseOutside = -1` sets marker phase to -1 | Box type uses PhaseInside/PhaseOutside incl. -1 assigned to PH2 | Seems indeed a bug; Please fix |
| S-021 | src/phase_transition.cpp:85, 100-119, 227-259, 426-433; dike.cpp:163 | Typos in phase transition input are silently accepted | extended Type/Parameter chains, still without final else |  |
| S-025 | src/adjoint.cpp:2348-2355, 2393-2400, 2439-2446 (llproX/Y/Z '=' assignments; xini copies 2370-2377, 2415-2422, 2461-2468) | Observation points in the same cell overwrite each other's weights | routed: written by Georg Reuber |  |
| S-026 | src/adjoint.cpp:633-661 (if(ts) bounds, else Ubar/Lbar=par_val; Ab=1 at 640, 655); used at 1049-1058 | Bounds: 0 means "no bound", and one bound pins every other parameter | Introduced 'if(ts)' bound test with InitGuess defaults; 3a9e100d added Ab=1 trigger for BLMVM. |  |
| S-027 | src/adjoint.cpp:1159 (fcconv size maxit), 1438 (descent write), 1534 (TAO write), 1456/1539 (count checks), 1094-1099 (print) | `fcconv` is written past its end | Moved fcconv creation into AdjointVectorsCreate; size unchanged. |  |
| S-028 | src/objFunct.h:85 (factor2array[51]); uses at src/adjoint.cpp:1296, 1328, 1343, 1356, 1415, 1445 | `factor2array[51]` is indexed up to 100 parameters | routed: written by Georg Reuber |  |
| S-031 | src/constEq.cpp:1008-1060 (setDataPhaseDiagram; ind[] at 1047-1050, clamps 1019-1058) | Phase diagram lookup reads one pressure row too low; edges wrong | Rewrote variable access (nT/nP, minP/dP); kept Georg's index logic. |  |
| S-037 | src/phase_transition.cpp:1302-1304 (zbounds[0],[1] test), 1322-1323, 1334; read at 413 | NotInAirBox uses only the first segment's `zbounds` | Original single-box NotInAirBox with temperature profiles using one z range. |  |
| S-039 | src/options.cpp:937-960 (setSolverOptions clear -> defaults -> reinsert); defaults building scaling/grid at 120-165, 251-256 | Command-line options are ignored while solver defaults are set | Original clear-options / set-defaults / reinsert-command-line structure. |  |
| S-040 | src/adjoint.cpp:2357-2380 (xini = Ae/f, f=v/v_node; Vy/Vz at 2402-2426, 2448-2472), 1650-1676 | Point-velocity cost isn't the point misfit, and its adjoint gradient doesn't match it | routed: written by Georg Reuber |  |
| S-041 | src/adjoint.cpp:868-875 (Ae = ts unscaled), 3455-3582 (derivatives, E2nd 'not yet implemented'); src/JacRes.cpp:675-677 | Strain-rate targets aren't scaled, and the strain-rate derivative is incomplete | Strain-rate target stored unscaled; incomplete (non-deviatoric) derivatives; scaling added only to printed value. |  |
| S-042 | src/adjoint.cpp:1744-1749 (inside lrank == 13 block) | `Adjoint_AdvectPoint` moves points on one CPU only, and again in every forward run | Restructured block (current blame owner); behaviour unchanged. |  |
| S-047 | src/adjoint.cpp:2505-2676 (Ap == 3 block; K lookup 2548, 2605, 2661); FindPointInCellAdjoint 102-125 | Surface observations (`Ap = 3`) use wrong grid positions and crash with several CPUs in z | routed: written by Georg Reuber |  |
| S-048 | src/adjoint.cpp:1152-1156 (Lb, Ub, val, P, grad with _MAX_PAR_ local size); 1186-1189 (aop->vx..sty with mdI local size) | Parameter vectors have 100 entries per CPU | routed: written by Georg Reuber |  |
| M-002 | src/JacRes.cpp:565-600 (JacResGetPShift, lines 582,599) | Top pressure shift isn't area-weighted and overwrites `p_shift` | equal-weight top-layer pressure average; overwrites ctrl.pShift |  |
| M-006 | src/surf.cpp:1711-1754 (FreeSurfSetInitialPerturbation) | Initial surface perturbation options have unit and NaN problems | added -FreeSurf_Wavelength/AmplCos/AmplNoise options; only AmplCos scaled |  |
| M-016 | src/bc.cpp:675-679, 712-713, 1775-1802 | Typos in BC options are silently ignored; `exy` drops no-slip on vx | bvel_face strcmp chain; unknown string leaves face=0 |  |
| M-039 | src/paraViewOutPassiveTracers.cpp:32 | `out_ptr` overwrites the `Passive_Tracer` switch after tracers are set up | added out_ptr read directly into ctrl.Passive_Tracer |  |
| M-042 | src/adjoint.cpp:937-947 | `Adjoint_ScaleCostFunction = Var` uses only the last point | routed: written by Georg Reuber |  |
| M-043 | src/adjoint.cpp:3026,3043,3109 (ScalingLaws) | Scaling-law exponents are wrong for `log10` parameters | b = grad*P/F scaling law; log10 parameters added later without adapting it |  |
| M-044 | src/adjoint.cpp:1328,1866,1995 | A parameter value of 0 gives a zero FD step | relative FD steps Perturb = CurVal*eps (zero for zero value) |  |
| M-062 | src/marker.cpp:1014,1152,1215; src/phase_transition.cpp:1197,1339 | Halfspace temperatures use hard-coded kappa = 1e-6 | hard-coded kappa=1e-6 for halfspace T (layer, box, PT box) |  |
| M-064 | src/phase_transition.cpp:886-905,922,1024-1029 | `ResetParam = APS` is applied every step and inconsistently | APS reset every step on reset side; PhaseOutside=-1 branch |  |
| M-065 | src/phase_transition.cpp:851-864 (BoxVicinity), 398-412 (bounds) | `BoxVicinity` uses wrong or stale bounds | BoxVicinity with dzBox computed from y bounds |  |
| M-070 | src/adjoint.cpp:204-260,296-400 | `AllMaterialParameters` scan: writes past arrays, checks only the first block | scan loop without blockID++, unchecked ExcludePhase/array limit, par_str leak |  |
| M-071 | src/adjoint.cpp:1934-1935,1946-1947,1995-2024 | Adjoint solve leaves the forward KSP on the `as_` prefix and cells in a perturbed state | adjoint reuses SNES KSP with as_ prefix; perturbed FormResidual leaves cell state |  |
| M-072 | src/adjoint.cpp:1120-1132,62-64,2742; src/LaMEMLib.cpp | `SyntheticForwardRun` saves a zero vector; the `OFdef = 0` path doesn't exist | routed: written by Georg Reuber |  |
| N-003 | src/JacRes.cpp:92,150; src/phase.cpp:1499; src/constEq.cpp:375,795-801,891-894 | `eta_min` is not a cutoff | default eta_st = eta_min in phase.cpp global-override routine |  |
| N-006 | src/phase.cpp:1458 (rho_c twice), 1467 (plasticity print) | PrintMatProp prints wrong values under some labels | added PrintMatProp, rho_c printed twice; inserted eta_vp in format but args left eta_st,eta_vp (swapped) |  |
| N-012 | src/bc.h:263; src/bc.cpp:782,827,1683,2961 (tSPCList partly used in JacResTemp.cpp:415) | Wrong comments and messages, debug print, small clean-ups | bc.h comment Gaussian=0/Poiseuille=1 wrong vs code |  |
| N-015 | src/objFunct.cpp (all); src/adjoint.cpp:3684-4150 | Dead code: whole ObjFunct module and the pointwise FD kernels | routed: written by Tobias Baumann, Georg Reuber |  |
| N-017 | src/phase.cpp:840-1352 (SetDiffProfile/SetDislProfile) | Creep profile experiment types not checked | added plagioclase profiles with experiment types |  |
| N-019 | src/nlsolve.cpp:209-232 (switch logic), 247 (MMFD label) | Newton/Picard switch logic | removed refRes=nrm at switch; divergence test vs first-iteration residual |  |
| N-021 | src/nlsolve.cpp:79-93,286; src/multigrid.cpp:283-284,1603-1802; src/multigrid.h:88; src/matData.cpp:159; src/lsolve.cpp:82-90; src/matAux.cpp:65,124-127 | Option overrides, stale comments, dead code, small clean-ups | MatMFFDSetType(WP) overrides -fd_mat_mffd_type |  |
| N-022 | src/AVD.cpp:98,403-665; src/advect.h:71; src/subgrid.cpp:740-786; src/cvi.cpp:13,157; src/marker.cpp:974,1071-1082; advect.cpp:1137; cvi.cpp:683,696 | Dead code, wrong comments and log lines, small clean-ups | marker.cpp layer comment says sin, code uses cos |  |
| N-025 | src/tssolve.cpp:216,234-296; src/options.cpp:906-907; src/scaling.cpp:118,183,241; src/LaMEMLib.cpp:986-1026; doc/src/man/GettingStarted.md:49-60 | Output and time-step quirks, labels, docs | nstep_out<=0 disables output; docs say nstep_out sets number of steps | the obvious ones can be fixed; nstep<=0 should indeed switch off output |
| N-027 | src/adjoint.cpp:449-458,804,903,951,1609,1701,1768-1776,2184-2199,2879,4155-4187 | Adjoint input and small quirks | optional-but-required Gradient/Scale params; stack Scaling pointer in modified DB; Adjoint_ApplyBCs sets dF to BC values | try fixing and see if tests still work |
| N-028 | src/JacRes.cpp:1086-1103 | FSSA only stabilizes the diagonal terms | restored diagonal-only form after trying full v-sum | leave as-is, as this is the key contribution. I'm pretty sure I had an option to do the full terms @ some stage |

### Main introducer over all 150 findings (not asked, for information)

Anton Popov: 52, Boris Kaus: 18, Andrea Piccolo: 14, Georg Reuber: 11, Garrett Ito: 9, Adina Püsök: 8, Jana Schierjott: 8, Arne Spang: 7, Yun Luo: 6, Tobias Baumann: 4, Wenrong Cao: 4, Iskander Ibragimov: 3, Panupong Kongpet: 2, Xinxin Wang: 1, Arthur Bauville: 1, Thomas Morrow: 1, Jianfeng Yang ("Jeff"): 1

Authors outside the review group and not routed (Panupong Kongpet, Iskander Ibragimov, Jianfeng Yang, …) are listed under "Introduced by" but are not asked.

## How the blame was done — caveats

- The audit was run against commit `6124aa7d`, which is not in this repository. Line numbers in the audit have drifted, so each finding was re-located at HEAD (c862179f) by content; the "At HEAD" line gives the current location. All findings were found at HEAD; none are fixed (M-019 is partly mitigated by `-snes_max_funcs unlimited` in options.cpp; N-017 is a physics question blame cannot answer).
- `git blame -w -M -C` with an ignore list of pure reformatting commits (astyle, tab cleanups), and where blame still landed on a bulk reformat/move, the introducing commit was traced with `git log -S/-L`. For "missing check" bugs, the line blamed is the code that lacks the check. Some attributions are judgement calls (e.g. S-003, S-035, S-039, S-042, M-046, M-067); the notes say so.
- Names are merged from aliases: lapopov/sgeotech/Anton Popov; APiccolo89/piccolo/Andrea Piccolo; Garrett/GIto/Garrett Ito; wenrongcao/Wenrong Cao; Yun Luo is `yunluo132@gmail.com` (single commit e4778c61, 2026-08-02, which added the FastScape coupling). **Please check:** "tosaba <baumann@uni-mainz.de>" was read as Tobias Baumann; "Jeff <yangjf@mail.iggcas.ac.cn>" as Jianfeng Yang.
- Commit dates come from git as recorded (many 2026 dates).
- The messages to the reviewers have not been sent — this file is the request, ready to be forwarded.

---

# Findings

All checked issues from the review, in three groups by severity: S-001… (major), M-001… (minor),
N-001… (note). Within each group, issues are sorted by effort: quick, medium, design.
The format is described in `README.md`.

Next free IDs: S-049, M-074, N-030

## Summary

### Severe (major)

| ID | Severity | Confidence | Kind | Effort | Where | Title | Status | Introduced by |
|----|----------|------------|------|--------|-------|-------|--------|---|
| S-001 | major | confirmed | memory | quick | JacRes.cpp:1524 | `init_pres` without `fix_phase` reads `phRat[-1]` | fixed | Anton Popov |
| S-002 | major | possible | physics | quick | constEq.cpp:587-592 | Softening length scale `Lm` looks inverted (Le/Lm instead of Lm/Le) | open | Anton Popov |
| S-003 | major | confirmed | physics | quick | phase.cpp:1049 | Profile "Diabase-Caristan_1982" has activation volume Vn = 1 m^3/mol | fixed | Adina Püsök (+Popov) |
| S-004 | major | confirmed | logic | quick | surf.cpp:1444-1462 | Sediment model 2 uses values left over from the previous node | fixed | Tobias Baumann |
| S-005 | major | confirmed | restart | quick | fastscape.cpp:2910, 2747 | FastScape restart loses basement and silt fraction if their output is off | open | Yun Luo (+Kaus) |
| S-006 | major | confirmed | logic | quick | fdstag.cpp:550-551 | Grid stretching moves the global box about 0, not about `bg_ref_point` | fixed | Anton Popov |
| S-007 | major | confirmed | logic | quick | fdstag.cpp:727-737 | CFL time step ignores negative velocities on non-uniform grids | fixed | Anton Popov |
| S-008 | major | confirmed | physics | quick | bc.cpp:723, 1087, 1121 | `bvel_temperature_constant` is scaled with the time unit, not as a temperature | open | Andrea Piccolo (+Kaus) |
| S-009 | major | confirmed | logic | quick | bc.cpp:1680-1725, 2054 | User `bvel_velout` is always overwritten by the mass-balance value | fixed | Iskander Ibragimov (+Popov) |
| S-010 | major | confirmed | logic | quick | bc.cpp:692-697, 1690-1700 | `velin_num_periods = 1` sets the inflow velocity to zero and ignores net periods | fixed | Iskander Ibragimov (+Kongpet) |
| S-011 | major | confirmed | restart | quick | bc.cpp:1704-1715 | After restart, net-only inflow adds the net offset twice (function `static`) | fixed | Panupong Kongpet |
| S-012 | major | confirmed | physics | quick | bc.cpp:3125-3134 | 3D Gaussian plume inflow is not mass-conserving (missing radius² factor) | open | Boris Kaus |
| S-013 | major | confirmed | logic | quick | advect.cpp:1217-1228 | Garbage collection deletes the wrong markers when the AVD delete list is unsorted | fixed | Adina Püsök (+Popov) |
| S-014 | major | confirmed | parallel | quick | cvi.cpp:1305-1414 | STAGP interpolation treats processor boundaries as domain walls | open | Anton Popov |
| S-015 | major | confirmed | logic | quick | marker.cpp:2358 | Oblique ridge axis is wrong unless the front edge is at y = 0 | open | Jana Schierjott |
| S-016 | major | confirmed | input | quick | marker.cpp:984 | Layer `rand_ampl` is divided by the max phase ID instead of the length scale | open | Boris Kaus |
| S-017 | major | confirmed | input, memory | quick | marker.cpp:1996-2094 | Phase diagram load fails silently (density 0); size not checked against buffer | open | Georg Reuber (+Kaus) |
| S-018 | major | confirmed | input, memory | quick | tools.cpp:627; marker.cpp:1687 | Polygon stretch factors after the last control polygon are left unset | fixed | Arne Spang |
| S-019 | major | confirmed | logic | quick | advect.cpp:222-234; marker.cpp:210-260 | `rand_noiseGP` moves markers after cell mapping without remapping | fixed | Tobias Baumann |
| S-020 | major | confirmed | logic | quick | phase_transition.cpp:756-763,811-812,870 | Box with `PhaseOutside = -1` (as documented) sets marker phase to -1 and stops the run | open | Boris Kaus (+Piccolo) |
| S-021 | major | confirmed | input | quick | phase_transition.cpp:85,100-119,227-259,667-711 | Typos in `Type`, `Parameter_transition`, `Name_Clapeyron` and bad `ID` are silently accepted | open | Andrea Piccolo (+Kaus) |
| S-022 | major | confirmed | logic, memory | quick | dike.cpp:1312,1331 | Dike stress time averaging loops `i < sx+ny` instead of `sx+nx` | open | Garrett Ito |
| S-023 | major | confirmed | physics | quick | dike.cpp:463-479 | Dike heating scales the whole cell conductivity by the dike fraction | fixed | Garrett Ito (+Schierjott) |
| S-024 | major | confirmed | physics | quick | passive_tracer.cpp:560-561 | Passive tracer pressure and temperature are read half a cell too high | fixed | Andrea Piccolo (+Popov) |
| S-025 | major | confirmed | logic | quick | adjoint.cpp:2342-2470 | Observation points in the same cell overwrite each other's weights | open | Georg Reuber |
| S-026 | major | confirmed | logic | quick | adjoint.cpp:632-662 | Bounds: 0 means "no bound", and one bound pins every other parameter | open | Boris Kaus |
| S-027 | major | confirmed | memory | quick | adjoint.cpp:1159,1438,1534 | `fcconv` is written past its end (TAO not limited by `Inversion_maxit`; off by one without TAO) | open | Georg Reuber (+Kaus) |
| S-028 | major | confirmed | memory | quick | objFunct.h:85; adjoint.cpp:1294-1450 | `factor2array[51]` is indexed up to 100 parameters | open | Georg Reuber |
| S-029 | major | confirmed | physics | medium | JacResAux.cpp:244-321 | Lithostatic and pore pressure are zero in boundary ghost cells, so boundary edges get half the value | fixed | Anton Popov (+("Jeff")) |
| S-030 | major | possible | physics | medium | JacResTemp.cpp:520-522 | Adiabatic heating sees a fake pressure jump at boundaries with inflow or outflow | open | Andrea Piccolo |
| S-031 | major | confirmed | logic | medium | constEq.cpp:1008-1060 | Phase diagram lookup reads one pressure row too low, and its edges are wrong | open | Georg Reuber (+Kaus) |
| S-032 | major | confirmed | physics | medium | surf.cpp:1229-1253 | Slope-dependent erosion never erodes peaks and ridges (centered gradient) | open | Wenrong Cao |
| S-033 | major | possible | physics | medium | surf.cpp:1203-1211 | Slope-erosion sub-step limit isn't a stability limit | open | Wenrong Cao |
| S-034 | major | confirmed | logic | medium | LaMEMLib.cpp:765; fastscape.cpp:1599 | With FastScape, max-angle smoothing is undone every step | open | Yun Luo |
| S-035 | major | confirmed | logic | medium | AVD.cpp:1192-1199 | AVD deletion ignores the Voronoi volumes (`tclaimed` is always 0) | fixed | Anton Popov (+Püsök) |
| S-036 | major | confirmed | physics | medium | phase_transition.cpp:592-594,667-711 | Clapeyron laws are wrong in SI units (C read as K, MPa/K read as Pa/K) | open | Arne Spang (+Piccolo) |
| S-037 | major | confirmed | logic | medium | phase_transition.cpp:1303,1322,1334 | NotInAirBox uses only the first segment's `zbounds` | open | Garrett Ito (+Kaus) |
| S-038 | major | confirmed | logic | medium | passive_tracer.cpp:585-607, 904-938, 1032-1049 | Passive tracers use a stale marker-to-cell list (after marker exchange) | open | Andrea Piccolo |
| S-039 | major | confirmed | input | medium | options.cpp:944-950, 126-162, 251-252 | Command-line options are ignored while solver defaults (MG levels, grid, thermal) are set | fixed | Anton Popov (+Kaus) |
| S-040 | major | confirmed | physics | medium | adjoint.cpp:2357-2380,1639-1663 | Point-velocity cost isn't the point misfit, and its adjoint gradient doesn't match it | open | Georg Reuber |
| S-041 | major | confirmed | physics | medium | adjoint.cpp:860-880,3455-3520 | Strain-rate targets aren't scaled, and the strain-rate derivative is incomplete | open | Boris Kaus (+Reuber) |
| S-042 | major | confirmed | parallel | medium | adjoint.cpp:1744-1749 | `Adjoint_AdvectPoint` moves points on one CPU only, and again in every forward run | open | Georg Reuber (+Kaus) |
| S-043 | major | possible | physics | design | fastscape.cpp:458-469, 646 | FastScape grid doesn't follow a stretching LaMEM box | open | Yun Luo |
| S-044 | major | confirmed | parallel | design | bc.cpp:2595-2604, 2648-2657 | `fix_phase` / `fix_cell` lose faces on processor boundaries; results depend on CPU count | open | Anton Popov |
| S-045 | major | confirmed | logic | design | advect.cpp:278,1086-1096 | With `exy_num_periods`, outflow markers are deleted, not wrapped | open | Anton Popov |
| S-046 | major | possible | physics | design | LaMEMLib.cpp:715; tssolve.cpp:308 | Solve and marker advection use different dt when dt is not fixed | open | Anton Popov |
| S-047 | major | confirmed | parallel, memory | design | adjoint.cpp:2505-2690,102-112 | Surface observations (`Ap = 3`) use wrong grid positions and crash with several CPUs in z | open | Georg Reuber |
| S-048 | major | possible | parallel | design | adjoint.cpp:1152-1157,1186-1189 | Parameter vectors have 100 entries per CPU, so TAO sees a vector that grows with CPU count | open | Georg Reuber |

### Minor

| ID | Severity | Confidence | Kind | Effort | Where | Title | Status | Introduced by |
|----|----------|------------|------|--------|-------|-------|--------|---|
| M-001 | minor | confirmed | logic | quick | JacResTemp.cpp:525,534 | Adiabatic heating efficiency is applied twice | open | Andrea Piccolo |
| M-002 | minor | confirmed | physics | quick | JacRes.cpp:558-602 | Top pressure shift isn't area-weighted and silently overwrites `p_shift` | open | Boris Kaus |
| M-003 | minor | confirmed | logic | quick | constEq.cpp:734-742 | `Phase_Melt` density ignores the `mfmax` cap | open | Andrea Piccolo |
| M-004 | minor | confirmed | logic | quick | phase.cpp:223-232 | Healing defaults are set only when printing is on | open | Garrett Ito |
| M-005 | minor | confirmed | logic | quick | surf.cpp:1531 | Sediment model 3 assumes the box starts at x = 0 | fixed | Tobias Baumann |
| M-006 | minor | confirmed | input, parallel | quick | surf.cpp:1711-1754 | Initial surface perturbation: unscaled inputs, NaN with noise only | open | Boris Kaus |
| M-007 | minor | confirmed | physics | quick | surf.cpp:1651, 1662 | Topographic diffusion boundary is not zero-flux | open | Wenrong Cao |
| M-008 | minor | confirmed | logic | quick | surf.cpp:1004, 1068, 1346 | Erosion and sedimentation levels overshoot by up to one step | fixed | Xinxin Wang (+Cao, Baumann) |
| M-009 | minor | possible | logic | quick | surf.cpp:461 | Surface outside the box silently gets zero velocity | fixed | Anton Popov |
| M-010 | minor | confirmed | memory | quick | fastscape.cpp:1234-1317, 3100 | Memory leaks in FastScape code | open | Yun Luo |
| M-011 | minor | possible | logic | quick | fastscape.cpp:1055 | Non-uniform FastScape interpolation can skip a node | open | Yun Luo |
| M-012 | minor | confirmed | logic, memory | quick | tools.cpp:628 | `interpStretch` can leave stretch factors unset after the last control polygon | fixed | Arne Spang |
| M-013 | minor | confirmed | input | quick | bc.cpp:848 | Phase 0 can't be used as permeable / plume mantle inflow phase | open | Andrea Piccolo |
| M-014 | minor | confirmed | logic | quick | bc.cpp:2097-2099, 2161-2163 | `relax_dist = 0` with `bvel_face_out` gives NaN when a cell center is exactly at top/bot | open | Andrea Piccolo |
| M-015 | minor | confirmed | physics | quick | bc.cpp:1480-1519, 3005-3017 | Plume temperature: different Gaussian width in BC and on injected markers; unset Tbot used as -1 | open | Panupong Kongpet (+Piccolo) |
| M-016 | minor | confirmed | input | quick | bc.cpp:675-679, 712-713, 1784-1802 | Typos in BC options are silently ignored; `exy` silently drops no-slip on vx | fixed | Anton Popov (+Kaus, Piccolo) |
| M-017 | minor | confirmed | logic | quick | matData.cpp:497-498, 529-530, 561-562 | Coarse edge viscosity at upper boundaries averages in zero ghost cells | fixed | Anton Popov |
| M-018 | minor | confirmed | logic | quick | nlsolve.cpp:430-515 | Some SNES/KSP failure reasons are not printed or not treated as severe | fixed | Anton Popov |
| M-019 | minor | possible | logic | quick | nlsolve.cpp:284 | MFFD mat-vecs count as SNES function calls; long Newton runs may hit max_funcs | won't fix | Anton Popov |
| M-020 | minor | confirmed | memory | quick | advect.cpp:1485-1539 | ADVCheckCorners reads an empty buffer when the cell and its local neighbours are empty | fixed | Adina Püsök |
| M-021 | minor | confirmed | logic | quick | AVD.cpp:753-754,1136-1137 | AVD never refills empty volumes and at most doubles sparse ones | not an issue | Anton Popov |
| M-022 | minor | confirmed | logic | quick | AVD.cpp:830-833 | `FindPointInCell` returns an error code as a cell index | fixed | Anton Popov |
| M-023 | minor | confirmed | logic | quick | cvi.cpp:204-205 | RK2 deletes markers whose midpoint leaves the domain | fixed | Adina Püsök |
| M-024 | minor | confirmed | logic | quick | subgrid.cpp:673; LaMEMLib.cpp:731-734 | `ADVGetSedPhase` uses stale cell numbers when called for passive tracers | open | Andrea Piccolo |
| M-025 | minor | confirmed | memory | quick | marker.cpp:604-642 | Temperature-file setup loops over the wrong marker count; index can run out of range | fixed | Arthur Bauville (+Popov) |
| M-026 | minor | confirmed | memory | quick | marker.cpp:294 | `ADVMarkSave` crashes when the save path has no '/'; APS not saved | partly fixed | Anton Popov |
| M-027 | minor | confirmed | logic | quick | advect.cpp:1089 | Periodic x-shift uses `>` while the rank test uses `>=` | fixed | Anton Popov |
| M-028 | minor | confirmed | physics | quick | advect.cpp:2053; bc.cpp:2952 | Adiabatic gradient uses \|z - Ztop\|, so rock above the reference level is heated | open | Andrea Piccolo |
| M-029 | minor | confirmed | input | quick | marker.cpp:1193-1200 | Ridge spreading rate can silently be 0 | open | Jana Schierjott |
| M-030 | minor | confirmed | input | quick | marker.cpp:378; advect.cpp:156-180; cvi.cpp:1239 | Marker input not checked (negative phase, nmark_sub, nmin > nmax, avd <= 0, stagp_a = 0) | fixed | Anton Popov (+Baumann, Püsök) |
| M-031 | minor | confirmed | memory | quick | cvi.cpp:800,826-828 | Memory leak of `vi->markind` in `ADVelReAllocStorage` | fixed | Adina Püsök |
| M-032 | minor | confirmed | memory | quick | dike.cpp:352,673 | Dikes without free surface read `phRat[-1]` | open | Garrett Ito |
| M-033 | minor | possible | logic | quick | dike.cpp:753-754,1095 | Dynamic dike can produce NaN (no neighbour within `filty`, or zero lithosphere thickness) | open | Garrett Ito |
| M-034 | minor | confirmed | restart | quick | dike.cpp:219,1588,1593 | Restart with dikes re-reads time stepping and resets the averaging counter | open | Garrett Ito |
| M-035 | minor | confirmed | input | quick | dike.cpp:212,306-345,1230 | Dike opening silently 0 without `bvel_velin`; `magPwidth` not scaled | open | Jana Schierjott (+Ito) |
| M-036 | minor | confirmed | logic, memory | quick | passive_tracer.cpp:599, 370, 607, 1049 | Tracer melt phase search checks the wrong marker; empty search reads past the list | open | Andrea Piccolo (+Luo) |
| M-037 | minor | confirmed | input | quick | tssolve.cpp:56, 61 | CFL and CFLMAX range checks can never trigger | fixed | Anton Popov |
| M-038 | minor | confirmed | parallel | quick | options.cpp:482 | Automatic coarse-grid CPU count uses the local, not the total, cell count | fixed | Anton Popov |
| M-039 | minor | confirmed | logic | quick | paraViewOutPassiveTracers.cpp:32 | `out_ptr` overwrites the `Passive_Tracer` switch after tracers are set up | open | Boris Kaus (+Piccolo) |
| M-040 | minor | confirmed | memory | quick | paraViewOutBin.cpp:770-772 | `UpdatePVDFile` uses the file name after `free` in its error message | fixed | Anton Popov |
| M-041 | minor | possible | restart | quick | paraViewOutBin.cpp:766-797 | After a restart from an older checkpoint, .pvd files keep old entries after `</VTKFile>` | fixed | Anton Popov |
| M-042 | minor | confirmed | logic | quick | adjoint.cpp:920-947 | `Adjoint_ScaleCostFunction = Var` uses only the last point | open | Georg Reuber |
| M-043 | minor | confirmed | physics | quick | adjoint.cpp:3024-3040,3080-3095 | Scaling-law exponents are wrong for `log10` parameters | open | Boris Kaus |
| M-044 | minor | confirmed | logic | quick | adjoint.cpp:1866,1995,1328 | A parameter value of 0 gives a zero FD step (division by zero) | open | Boris Kaus (+Reuber) |
| M-045 | minor | confirmed | parallel | medium | JacRes.cpp:1004-1017 | Dike strain-rate correction is missing in ghost cells, so results depend on CPU count | open | Jana Schierjott |
| M-046 | minor | possible | physics | medium | constEq.cpp:798-822 | Dike opening still enters the stabilization stress and shear heating | open | Jana Schierjott (+Popov) |
| M-047 | minor | possible | logic | medium | JacRes.cpp:1545-1703 | `init_lith_pres` only sets the pressure history, and its stop test is weak | open | Anton Popov (+Spang) |
| M-048 | minor | confirmed | input | medium | phase.cpp:798-803; marker.cpp:523 | Phase temperature `T = 0` (°C) is treated as "not set" | open | Anton Popov (+Baumann, Ito) |
| M-049 | minor | possible | physics | medium | constEq.cpp:603-626, 697, 756 | Shear and bulk moduli averaged with non-elastic phases counted as 0 | open | Anton Popov |
| M-050 | minor | confirmed | logic | medium | surf.cpp:1167, 1223, 1643 | Slope erosion and topographic diffusion ignore periodic x | open | Wenrong Cao |
| M-051 | minor | confirmed | parallel, restart | medium | tools.cpp:80, 293-322; fdstag.cpp:292-294 | Rank-0-only errors can hang; restart reads aren't checked | fixed | Anton Popov |
| M-052 | minor | confirmed | input | medium | bc.cpp:2103-2123, 2167-2187 | Sign of `velin` on Right/Back faces flips depending on `bvel_face_out` | open | Iskander Ibragimov (+Piccolo) |
| M-053 | minor | confirmed | parallel, input | medium | bc.cpp:1229-1274, 1148 | Fixed-cell file read: unchecked fread/stat, one-rank errors hang, file not closed on error | fixed | Anton Popov |
| M-054 | minor | possible | input | medium | bc.cpp:2132-2142, 2213-2223 | CompensatingInflow has no mass-balance check (needs user `velbot`/`veltop`) | fixed | Thomas Morrow (+Kongpet) |
| M-055 | minor | confirmed | logic | medium | matFree.cpp:637-684, 1304-1309, 1355-1360; matBFBT.cpp:208-209 | Matrix-free operator and diagonal treat the periodic x seam as a free-slip wall | fixed | Anton Popov |
| M-056 | minor | confirmed | logic | medium | matData.cpp:642-645, 808-811 | Rediscretized coarse levels lose no-slip and pressure BCs on boundary ghosts | partly fixed | Anton Popov |
| M-057 | minor | confirmed | parallel | medium | advect.cpp:1368-1568 | ADVCheckCorners depends on CPU count, reuses the same seed, misses markers on the centre plane | fixed | Adina Püsök (+Popov) |
| M-058 | minor | confirmed | parallel | medium | AVD.cpp:747-749,791-793 | AVD half-volume rule uses processor boundaries; nmax not halved | fixed | Anton Popov |
| M-059 | minor | confirmed | logic | medium | AVD.cpp:813-816,865-884,944-946 | AVD "edge" control volumes are face volumes; printed labels swapped | fixed | Anton Popov |
| M-060 | minor | confirmed | memory | medium | marker.cpp:1466,1622-1744 | Polygon file sizes not checked against buffers | fixed | Arne Spang |
| M-061 | minor | confirmed | logic | medium | interpolate.cpp:102,158,220,377 | Corner output at the periodic x seam uses one side only | fixed | Anton Popov |
| M-062 | minor | confirmed | physics | medium | marker.cpp:1013,1152,1215; phase_transition.cpp:1197,1339 | Halfspace temperatures use hard-coded kappa = 1e-6 (wrong in non-dim mode) | open | Boris Kaus (+Schierjott) |
| M-063 | minor | confirmed | parallel | medium | subgrid.cpp:85,606-608; advect.cpp:1906; tools.cpp:293,798 | Rank-local hard errors can hang; sedimentation-phase error on fast uplift | fixed | Anton Popov (+Püsök) |
| M-064 | minor | confirmed | logic | medium | phase_transition.cpp:883-905,922,1026-1029 | `ResetParam = APS` resets every step; side depends on type; broken for NotInAirBox | open | Boris Kaus (+Schierjott) |
| M-065 | minor | confirmed | logic | medium | phase_transition.cpp:851-864 | `BoxVicinity` uses y-size for z, and stale or zero bounds for NotInAirBox | open | Boris Kaus (+Ito) |
| M-066 | minor | confirmed | memory | medium | tssolve.cpp:381-423, 457-462 | Time-step schedule can overflow its 2000-entry arrays and cuts nstep_max | fixed | Arne Spang |
| M-067 | minor | confirmed | logic | medium | LaMEMLib.cpp:713, phase_transition.cpp:967 | Time-step restart (CFLMAX) repeats phase transitions, so moving boxes move twice | open | Jana Schierjott (+Ito, Popov) |
| M-068 | minor | confirmed | logic, input | medium | passive_tracer.cpp:262, 420-797, 63-91 | Passive tracer APS never updated, IDs wrong when nx != ny, typos silently mean "Always" | open | Andrea Piccolo (+McMillan) |
| M-069 | minor | confirmed | physics | medium | outFunct.cpp:742-743, 767-768; JacResAux.cpp | SHmax/EHmax direction vectors are averaged to corners, which breaks near north-south | fixed | Anton Popov |
| M-070 | minor | confirmed | memory | medium | adjoint.cpp:204-389 | `AllMaterialParameters` scan: writes past arrays, checks only the first block | open | Boris Kaus |
| M-071 | minor | confirmed | logic | medium | adjoint.cpp:1935,1947,1995-2024 | Adjoint solve leaves the forward KSP on the `as_` prefix and cells in a perturbed state | open | Boris Kaus (+Reuber) |
| M-072 | minor | confirmed | cleanup | medium | adjoint.cpp:1123-1132,62-66 | `SyntheticForwardRun` saves a zero vector; the `OFdef = 0` path doesn't exist | open | Georg Reuber |
| M-073 | minor | confirmed | parallel | design | dike.cpp:111-120,959-1001,1314,1395 | Dynamic dike parallel limits are not checked up front; results depend on CPU layout | open | Garrett Ito |

### Notes

| ID | Severity | Confidence | Kind | Effort | Where | Title | Status | Introduced by |
|----|----------|------------|------|--------|-------|-------|--------|---|
| N-001 | note | confirmed | logic | quick | JacRes.cpp:917 | `dikeRHS` is not set when dikes are off | open | Jana Schierjott |
| N-002 | note | confirmed | cleanup | quick | JacResAux.cpp:389; JacResTemp.cpp:85,107 | Material data changed inside compute loops; Controls copied for every cell | fixed | Jana Schierjott (+Popov) |
| N-003 | note | confirmed | physics | quick | phase.cpp:1511; constEq.cpp:375,798 | `eta_min` is not a cutoff; it's added to every phase's viscosity | open | Anton Popov (+Kaus) |
| N-004 | note | confirmed | input | quick | phase.cpp:526-548 | Creep input checks never fire; `Bd`/`Bn` are silently overwritten by `eta`/`eta0` | open | Anton Popov (+Kongpet) |
| N-005 | note | confirmed | cleanup | quick | meltParam.cpp | Whole file is dead code (no callers) | open | Jianfeng Yang ("Jeff") (+Piccolo) |
| N-006 | note | confirmed | logic | quick | phase.cpp:1471,1479 | PrintMatProp prints wrong values under some labels | open | Boris Kaus |
| N-007 | note | confirmed | cleanup | quick | surf.cpp, fastscape.cpp | Small clean-ups (duplicate reads, wrong comments/units, folder name rounding) | fixed | Yun Luo (+Cao, Baumann) |
| N-008 | note | confirmed | physics | quick | surf.cpp:288-292 | Average topography is a plain node average | fixed | Anton Popov |
| N-009 | note | confirmed | cleanup | quick | Tensor.cpp:112-318, 493-603 | ISA / eigenvalue code is unused and has latent bugs | fixed | Anton Popov |
| N-010 | note | confirmed | cleanup | quick | fdstag.cpp:80-83 | Non-unit bias is rejected, so all bias code is dead | fixed | Anton Popov (+Luo) |
| N-011 | note | confirmed | input | quick | bc.cpp:2405-2411, 2234, 2369 | Velocity cylinder/box edge cases: NaN axis, skipped at initial guess unlike other BCs | fixed | Arne Spang (+Piccolo) |
| N-012 | note | confirmed | cleanup | quick | bc.h, bc.cpp (various) | Wrong comments and messages, debug print every step, small clean-ups | open | Andrea Piccolo (+Kaus, Kongpet) |
| N-013 | note | possible | logic | quick | AVD.cpp:1123-1125 | AVD point with no cells gives NaN coordinates | fixed | Adina Püsök (+Popov) |
| N-014 | note | confirmed | input | quick | input_file.dat:747-749; phase_transition.cpp:562-580 | Docs say Clapeyron parameters overwrite a predefined law; they are ignored | fixed | Anton Popov (+Spang) |
| N-015 | note | confirmed | cleanup | quick | objFunct.cpp; adjoint.cpp:3684-4150 | Dead code: whole ObjFunct module and the pointwise FD kernels | fixed | Tobias Baumann (+Reuber) |
| N-016 | note | confirmed | physics | medium | JacResAux.cpp:357 | `gw_level_type = surf` gives a flat water table at the mean topography | open | Anton Popov |
| N-017 | note | possible | physics | medium | phase.cpp:888-1352 | Experiment type of the creep profiles not checked against the sources | not an issue | Anton Popov (+Kaus) |
| N-018 | note | confirmed | memory | medium | tools.cpp:608-609, 658-662 | Polygon helpers use stack arrays and can divide by zero | fixed | Arne Spang |
| N-019 | note | confirmed | logic | medium | nlsolve.cpp:216-238 | Newton/Picard switch: divergence test uses initial residual; max-Newton gives one Picard step | open | Boris Kaus (+Popov, Spang) |
| N-020 | note | confirmed | physics | medium | lsolve.cpp, matrix.cpp, matAux.cpp | `jp_pgamma` means different things per solver; inv_eta Schur with penalty drops the eta term | not an issue | Anton Popov |
| N-021 | note | confirmed | cleanup | medium | nlsolve, lsolve, multigrid, matData, matFree (various) | Option overrides, stale comments, dead test code, small clean-ups | partly fixed | Anton Popov (+Kaus) |
| N-022 | note | confirmed | cleanup | medium | marker, advect, AVD, cvi, subgrid (various) | Dead code, wrong comments and log lines, small clean-ups | partly fixed | Adina Püsök (+Popov, Kaus) |
| N-023 | note | confirmed | cleanup | medium | phase_transition.cpp, dike.cpp (various) | Dead code, misleading names and comments, small clean-ups | open | Garrett Ito (+Piccolo, Schierjott) |
| N-024 | note | confirmed | input | medium | parsing.cpp (various) | Input parsing: silent bad numbers, octal integers, dropped values, rank-0 error hang | fixed | Anton Popov |
| N-025 | note | confirmed | cleanup | medium | options.cpp, tssolve.cpp, scaling.cpp, LaMEMLib.cpp, docs | Output and time-step quirks, labels, docs | partly fixed | Anton Popov (+Kaus, Spang) |
| N-026 | note | confirmed | cleanup | medium | outFunct.cpp, paraViewOutBin.cpp, paraViewOutAVD.cpp, paraViewOutPassiveTracers.cpp | Output quirks: names, stale residual, dead vectors, unchecked inputs | partly fixed | Anton Popov (+("Jeff"), Thielmann) |
| N-027 | note | confirmed | input | medium | adjoint.cpp (several) | Adjoint input and small quirks | open | Boris Kaus (+Reuber) |
| N-028 | note | confirmed | physics | design | JacRes.cpp:1086-1103 | FSSA only stabilizes the diagonal terms | won't fix | Anton Popov (+Kaus) |
| N-029 | note | confirmed | physics | design | subgrid.cpp:528-606 | Rock/air phase change keeps old stress, APS and T | won't fix | Anton Popov |

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
- Introduced by (git blame at HEAD):
  - Anton Popov — 15a0094f Linear pressure initial guess between top and bottom boundaries (2017-07-26) — wrote phRat[fixPhase] test without checking fixPhase == -1
  - Anton Popov — 8679b765 matrix constraint (started) (2017-05-11) — set fixPhase default -1 in bc.c
- At HEAD: still present, src/JacRes.cpp:1524. Note: init_pres code (2017) relies on fix_phase being set; no cross-check added.
- Feedback requested from: Anton Popov
- Status: fixed (group A)
- Boris: please fix

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
- Introduced by (git blame at HEAD):
  - Anton Popov — 49845492 Stabilizations and regularizations. (2020-05-10) — introduced Lm-based APS1/APS2 *= Le/Lm scaling
- At HEAD: still present, src/constEq.cpp:588-592. Note: Single commit introduced it; audit rates it only 'possible'.
- Feedback requested from: Anton Popov
- Status: open
- Boris: Anton, please check (hard to believe that this is wring but could be)
 
### S-003 [major, confirmed] Profile "Diabase-Caristan_1982" has Vn = 1 m^3/mol
- Tags: major, confirmed, physics, quick
- Where: src/phase.cpp:1049
- What: activation volume `m->Vn = 1;`, which is m^3/mol. Normal values are around 1e-5.
  At 1 GPa, P*V/RT is about 1e5, so `exp(-Q)` is 0 and dislocation creep is switched off
  everywhere except at the surface. All the other profiles use 0 or ~1e-5.
- Fix idea: probably `0` (the source in the comment gives no V). Check the reference.
- Introduced by (git blame at HEAD):
  - Adina Püsök — c238857d added diffusion and dislocation creep profiles from the literature (2015-03-06) — originally added the profile with Vn = 1 (comment: varies 2-10e-6)
  - Anton Popov — 64d9ea3b . (2017-02-16) — rewrote/moved profiles into phase.c carrying Vn = 1 (blame target)
- At HEAD: still present, src/phase.cpp:1037. Note: Original value from Adina Puesoek in 2015 (matProps.c); Anton's 2017 file reorganisation is what blame shows.
- Feedback requested from: Anton Popov
- Status: fixed (group A)
- Boris: should be zero (thats how it is in GeoParams - see https://github.com/JuliaGeodynamics/GeoParams.jl/blob/32dfaedb509ee1258551bf59c12a159e250fdfa1/src/CreepLaw/Data/DislocationCreep.jl#L190)

### S-004 [major, confirmed] Sediment model 2 uses values left over from the previous node
- Tags: major, confirmed, logic, quick
- Where: src/surf.cpp:1444-1462
- What: `zprop` and `zpropn` are set only in some branches. If `dn<0`, only `zprop` is set and
  `zpropn` keeps the previous node's value. If `d<0`, it's the other way round. Then
  `dz = zpropn - zprop` mixes two nodes. At time 0, `z = zprop` can also be stale.
- Why it matters: wrong sediment thickness next to the margin, and the result depends on loop
  order, so also on how the domain is split across CPUs.
- Fix idea: set both values in every branch (probably both = InitLevel + hUp in the first two).
- Introduced by (git blame at HEAD):
  - Tobias Baumann — d7e274d2 changed sedimentation surface into sedimentation increment (2019-03-22) — added zprop/zpropn dn<0 / d<0 branches setting only one value
  - Tobias Baumann — 283d15db added another sedimentation process (2019-03-18) — created the sediment model 2 branch structure and variable declarations
- At HEAD: still present, src/surf.cpp:1444-1462. Note: Raw name 'tosaba' (baumann@uni-mainz.de) assumed to be Tobias Baumann.
- Feedback requested from: Anton Popov
- Status: fixed (group A)
- Boris: yes could be wrong, please fix

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
- Introduced by (git blame at HEAD):
  - Yun Luo — e4778c61 New coupling module based on LaMEM 3.0.1 (2026-08-02) — copy of basement/silt fraction conditional on output flag; restart always uses them
  - Boris Kaus — 03205112 Merge branch 'master' into yl/fastscape-coupling-new (2026-08-13) — merge resolution turned output copy into out_flags[] array loop (same behaviour)
- At HEAD: still present, src/fastscape.cpp:2776-2804 and 2905-2914. Note: Behaviour from Yun's original manage_output(out_flag,...); Boris's merge only restructured it. Defaults are out=1, so bug needs output turned off.
- Feedback requested from: Yun Luo, Boris Kaus
- Status: open
- Boris: dont think this was ever tested wiuth restarting; would be good to have 

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
- Introduced by (git blame at HEAD):
  - Anton Popov — 157753dd Bug fixes. (2017-03-24) — added gcrdbeg/gcrdend *= (1+eps) in Discret1DStretch (ref=0 assumed)
- At HEAD: still present, src/fdstag.cpp:550-551. Note: Node stretch with ref point from 28f94019 (2014) was later; global bounds ignore ref.
- Feedback requested from: Anton Popov
- Status: fixed (group A); test with non-zero `bg_ref_point` not yet added
- Boris: Please fix & add a test

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
- Introduced by (git blame at HEAD):
  - Anton Popov — 65fb1ee7 time stepping controls (belongs to previous update) (2014-11-03) — wrote idt = v/h for non-uniform grid (no abs value)
- At HEAD: still present, src/fdstag.cpp:726-737. Note: Boris's ad39e5e0 (2014-11-04 bugfixes) only touched the following 'if(idt > idtmax)' line.
- Feedback requested from: Anton Popov
- Status: fixed (group A)
- Boris: fix

### S-008 [major, confirmed] `bvel_temperature_constant` is scaled with the time unit, not as a temperature
- Tags: major, confirmed, physics, quick
- Where: src/bc.cpp:723 (read), 1087 (print), 1121 (nondimensionalize) (commit 6124aa7d)
- What: the value is read with `scal->time` as scale, so it is divided by the time unit. The
  print multiplies it back, so the log shows the right number. Then line 1121 treats the
  divided value as °C: `(T/time_unit + Tshift)/T_unit`. Other temperatures use scale 1.0.
- Why it matters: with `bvel_temperature_inflow = Constant_T_inflow` in GEO or SI units the
  inflow markers get a wrong temperature (close to 0 °C in GEO units).
- Fix idea: read with scale 1.0 and print without `*scal->time`.
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 9c06b04e Merge branch 'master' into andrea_piccolo/inflow_bug_fix (2020-11-12) — first appearance of bvel_constant_temperature read with scal->time scale
  - Boris Kaus — ad72a8d5 renaming bvel_face to used keywords... proper formatting of print statements (2020-11-16) — added print multiplying by scal->time, masking the wrong scale in the log
- At HEAD: still present, src/bc.cpp:723 (read), 1087 (print), 1123 (nondimensionalize). Note: Read line first appears in Andrea's merge commit (not in either parent); renamed keyword in 073d032c (Andrea). Pickaxe skips merges.
- Feedback requested from: Andrea Piccolo, Boris Kaus
- Status: open
- Boris: hard to believe that this is wrong..

### S-009 [major, confirmed] User `bvel_velout` is always overwritten by the mass-balance value
- Tags: major, confirmed, logic, quick
- Where: src/bc.cpp:1680, 1699, 1719, 1725 (BCGetVelins), 2054 (BCApplyBoundVel)
- What: BCApplyBoundVel calls BCGetVelins first, and every branch of BCGetVelins sets
  `velout = -velin*(top-bot)/(bot-bz)`. The `bvel_velout` input read in BCCreate is never used.
- Why it matters: a documented input has no effect, silently.
- Fix idea: keep a flag for "user velout given" and only recompute when it is not set.
- Introduced by (git blame at HEAD):
  - Iskander Ibragimov — d6e7958a Added new inflow from boundary condition options (2023-05-17) — added BCGetVelins call each step, recomputing velout in every branch
  - Anton Popov — 5d3bb5f1 boundary velocity constraint implementation (2015-11-07) — original mass-balance velout formula (previously only computed if velout unset)
- At HEAD: still present, src/bc.cpp:1679-1725 (BCGetVelins), 2047 (BCApplyBoundVel). Note: Before d6e7958a, velout computed in BCCreate only when bvel_velout was unset (DBL_MAX default).
- Feedback requested from: Anton Popov
- Status: fixed (group B, batch 1)
- Boris: I would anyways like to have a much more general way to set boundary & internal velocities. Adding a flag seems a good idea; default behavior should still be to compute this in a mass-conservative way

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
- Introduced by (git blame at HEAD):
  - Iskander Ibragimov — d6e7958a Added new inflow from boundary condition options (2023-05-17) — introduced 'if(VelNumPeriods)' branch while array read only for >1
  - Panupong Kongpet — 7d0eb7a1 Add velin_net_num_periods, velin_net_time_delims, bvel_velin_net (2026-01-25) — added net-only branch after the shadowing VelNumPeriods branch
- At HEAD: still present, src/bc.cpp:690-699, 1690-1700. Note: Merge b74d610e (Anton, 2026-04-20) only combined these; logic is theirs.
- Feedback requested from: Anton Popov (assigned by code area; none of the main suspects wrote these lines)
- Status: fixed (group A)
- Boris: please fix

### S-011 [major, confirmed] After restart, net-only inflow adds the net offset twice
- Tags: major, confirmed, restart, quick
- Where: src/bc.cpp:1704-1715
- What: with only `velin_net` periods, the base inflow is stored in a function-`static`
  the first time BCGetVelins runs, then `bc->velin = base + net[kk]` is written back. The
  whole BCCtx is saved in the restart file, so after restart `bc->velin` already holds
  `base + net`, and the fresh `static` takes that as the base. Result: `base + 2*net`.
- Why it matters: restarted runs use a different inflow velocity than the original run.
- Fix idea: store the base value in BCCtx (set once in BCCreate) instead of a `static`.
- Introduced by (git blame at HEAD):
  - Panupong Kongpet — 7d0eb7a1 Add velin_net_num_periods, velin_net_time_delims, bvel_velin_net (2026-01-25) — added function-static baseVelin captured from bc->velin
- At HEAD: still present, src/bc.cpp:1704-1718
- Feedback requested from: Anton Popov (assigned by code area; none of the main suspects wrote these lines)
- Status: fixed (group A)
- Boris: please fix if you think this makes sense

### S-012 [major, confirmed] 3D Gaussian plume inflow is not mass-conserving (missing radius² factor)
- Tags: major, confirmed, physics, quick
- Where: src/bc.cpp:3125-3134
- What: the integral of `exp(-r²/c²)` over the box is `(π c²/4)·Δerf_x·Δerf_y`. The code uses
  `π/4·Δerf_x·Δerf_y/Area`, without `c²`. The 2D branch has the `c` factor and is right.
- Why it matters: with `Plume_VelocityType = Gaussian` and `Plume_Dimension = 3D` the outflow
  velocity is wrong, so net volume flux through the bottom is not zero. In an
  incompressible closed box this has no solution or forces flow through other boundaries.
- Fix idea: multiply `a, b, d, e` by `c*c`.
- Introduced by (git blame at HEAD):
  - Boris Kaus — b9d87321 added a gaussian-type inflow/outflow velocity at the bottom (2020-10-20) — derived 3D erf integral without the c^2 factor
- At HEAD: still present, src/bc.cpp:3125-3134. Note: Arne Spang (df705292, 2022) only added areaFrac to V_out line.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: please fix; may explain some of the issue we had wuth this,

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
- Introduced by (git blame at HEAD):
  - Adina Püsök — 766f54c3 Marker control (solves 'Empty Control Volume' error) (2015-03-30) — AVD delete list idel filled by sorted volume/area order, not index
  - Anton Popov — 88b60a16 code integration & functionality restore update version 0.0 (2014-09-09) — garbage-collect loop that assumes ascending idel
- At HEAD: still present, src/advect.cpp:1217-1228; src/AVD.cpp:591,1205. Note: Loop correct for original ascending lists; bug arose when Adina added AVD marker deletion (2015).
- Feedback requested from: Anton Popov
- Status: fixed (group A)
- Boris: please fix even when his is no longer default in thye code

### S-014 [major, confirmed] STAGP interpolation treats processor boundaries as domain walls
- Tags: major, confirmed, parallel, quick
- Where: cvi.cpp:1305/1314 (`IN == 0`, `IN == nx`), 1355/1364 (`JN`), 1405/1414 (`KN`) (commit 6124aa7d)
- What: `nx = fs->dsx.ncels` and `IN` are local. At an internal processor boundary the code
  uses the wall extrapolation `2*v[IN]-v[IN+1]` instead of the available ghost values.
- Why it matters: `interp = stagp` gives different marker paths for different CPU counts.
- Fix idea: test global indices (`sx+IN == 0`, `sx+IN == global nx`). The comment "ONLY FOR
  FREE SLIP BC" (cvi.cpp:1298) is not enforced either.
- Introduced by (git blame at HEAD):
  - Anton Popov — 8763ccaf Multiple updates from concurrent branches (2017-01-11) — added STAGP interpolation with local IN==0/nx wall tests
- At HEAD: still present, src/cvi.cpp:1305/1314, 1355/1364, 1405/1414
- Feedback requested from: Anton Popov
- Status: open (moved from A to C: needs 2 ghost layers or another stencil)
- Boris: seems urgent to fix

### S-015 [major, confirmed] Oblique ridge axis is wrong unless the front edge is at y = 0
- Tags: major, confirmed, logic, quick
- Where: marker.cpp:2358 (commit 6124aa7d)
- What: `x_oblique = (xL-xR)/(yF-yB) * y + xL;` should use `(y - yF)`. Divides by zero if
  yF == yB.
- Why it matters: plate age and halfspace temperature are wrong in the whole ridge box.
  Tests only use straight ridges.
- Fix idea: use `(y - y_ridgeFront)` and check yF != yB.
- Introduced by (git blame at HEAD):
  - Jana Schierjott — 993280d8 option added for oblique and perpendicular ridges (2020-11-17) — wrote x_oblique without (y - yF) offset
- At HEAD: still present, src/marker.cpp:2358
- Feedback requested from: Garrett Ito
- Status: open
- Boris: yes, seems wronng; fix

### S-016 [major, confirmed] Layer `rand_ampl` is divided by the max phase ID
- Tags: major, confirmed, input, quick
- Where: marker.cpp:984, 975, 2207 (commit 6124aa7d)
- What: `getScalarParam(..., "rand_ampl", ..., 1, (PetscScalar) maxPhaseID)` divides by the
  last argument, which should be `chLen`. With one phase (maxPhaseID = 0) it gives inf.
  `cosine` is range-checked against maxPhaseID too, so `cosine = 2` is silently ignored.
  The noise uses `rand()` per marker, so it depends on CPU count.
- Why it matters: layer interface noise has the wrong size, or inf/NaN layer bounds.
- Fix idea: pass `chLen`; check `cosine` against 1; use PetscRandom.
- Introduced by (git blame at HEAD):
  - Boris Kaus — 01b07ee4 added a Rayleigh-Taylor Instability benchmark (2020-05-26) — added rand_ampl read with maxPhaseID instead of chLen; rand() noise
  - Boris Kaus — 8dd159de added cosine perturbation of layers (2020-05-22) — added cosine range-checked against maxPhaseID
- At HEAD: still present, src/marker.cpp:985, 976, 2207
- Feedback requested from: Boris Kaus
- Status: open
- Boris: not sure how this happened; fix please

### S-017 [major, confirmed] Phase diagram load fails silently; size not checked against buffer
- Tags: major, confirmed, input, memory, quick
- Where: marker.cpp:1996-2000, 2090-2094, 2060-2088 (commit 6124aa7d)
- What: more than 8 diagrams ("buffer too small") or an unknown column count only print a
  message and return. The name is not stored, so constEq.cpp:990 sets `pd->rho = 0`.
  `nT*nP` is never checked against `_max_pd_sz_` (40100), so a 201×201 file writes past
  `rho_v`. `fopen`/`fscanf` results are not checked; printed max T/P are one step too wide.
- Why it matters: runs continue with zero density, or memory is overwritten.
- Fix idea: SETERRQ in both cases; check n <= `_max_pd_sz_`; check file I/O.
- Introduced by (git blame at HEAD):
  - Georg Reuber — 61aa7159 Added phase diagrams to LaMEM (2018-01-19) — wrote loader with silent returns, no size/file checks
  - Boris Kaus — e92e969b corrected a bug in saving *.in file from Perple_X, more checks (2021-02-04) — wrote T/P range print (minor)
- At HEAD: still present, src/marker.cpp:1996-2000, 2088-2094, 2043 (and 2054-2086)
- Feedback requested from: Boris Kaus
- Status: open
- Boris: Please fix

### S-018 [major, confirmed] Polygon stretch factors after the last control polygon are left unset
- Tags: major, confirmed, input, memory, quick
- Where: tools.cpp:627 (`interpStretch`); marker.cpp:1687-1696 (commit 6124aa7d)
- What: the tail is filled only if `CtrlPoly[numCtrlPoly-1] != numCtrlPoly-1`; it should
  compare with `numPoly-1`. E.g. one control polygon at position 1: all later slices use
  uninitialised stack values (only clamped to >= 0.1). In marker.cpp, PolyPos < 1, wrong
  order and unknown VolID are not checked; PolyID is read but unused.
- Why it matters: polygon bodies get random sizes.
- Fix idea: always fill up to numPoly-1; check 1 <= Pos <= Vol.num, increasing, VolID valid.
- Introduced by (git blame at HEAD):
  - Arne Spang — a88df698 added geometry manipulation (2020-07-08) — wrote interpStretch tail fill using numCtrlPoly-1 instead of numPoly-1
  - Arne Spang — 82be9435 Extended ability to read geomIO polygons (2022-05-05) — polygon control reading without position/volume checks
- At HEAD: still present, src/tools.cpp:628-635; src/marker.cpp:1685-1697
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### S-019 [major, confirmed] `rand_noiseGP` moves markers after cell mapping without remapping
- Tags: major, confirmed, logic, quick
- Where: advect.cpp:222-234; marker.cpp:210-260 (commit 6124aa7d)
- What: `ADVMarkPerturb` runs after `ADVMapMarkToCells` and shifts markers by ±0.5 sub-cell.
  With `rand_noise` also on, markers can cross into the next cell or rank.
- Why it matters: `cellnum` is stale, so the first projection puts phase and T in the wrong
  cell; CheckMarkers misses it; markers off-rank cause "non-local marker" aborts.
- Fix idea: call `ADVMapMarkToCells` after the perturbation, or clamp to the cell.
- Introduced by (git blame at HEAD):
  - Tobias Baumann — b5152db0 added option to perturb markers after setting up geometric structure (2019-03-29) — added ADVMarkPerturb after ADVMapMarkToCells, no remapping
- At HEAD: still present, src/advect.cpp:225-228; src/marker.cpp:210-262. Note: Anton's d3a5c7ee (2019-04-02) only rearranged the surrounding calls.
- Feedback requested from: Anton Popov
- Status: fixed (group A): perturbed coordinate kept in host cell
- Boris: Please fix

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
- Introduced by (git blame at HEAD):
  - Boris Kaus — 105ac06b adding a Box-like phase transition (2020-11-16) — Box type uses PhaseInside/PhaseOutside incl. -1 assigned to PH2
  - Boris Kaus — 32272491 fix for phase transitions (2024-02-04) — nphc test line includes PhaseAbove/Below zeros for Box
  - Andrea Piccolo — d7aff650 Bug Fix: communicate phase transition not changing phase (2022-03-02) — introduced nphc flag
- At HEAD: still present, src/phase_transition.cpp:756-763, 811-812, 824-825. Note: Andrea's raw name 'andrea@Piccolo-W10.bgi.uni-bayreuth.de' also APiccolo89.
- Feedback requested from: Andrea Piccolo, Boris Kaus
- Status: open
- Boris: Seems indeed a bug; Please fix

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
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 672cf0c2 First attempt to do the phase transitions routine (2020-05-07); a97d96fc Return to a marker based approach (2020-07-03); 3cc44e46 change all the string comparison (2020-10-19) — ID range check and string-compare chains without else error
  - Boris Kaus — 4bc1f1a8 enum underscores, minor bugfixes for PhaseTransitions (2020-10-19); 570950e4 melt-fraction transition (2020-11-17) — extended Type/Parameter chains, still without final else
- At HEAD: still present, src/phase_transition.cpp:85, 100-119, 227-259, 426-433; dike.cpp:163. Note: Dike ID check (dike.cpp:163) not blamed (dike author not examined).
- Feedback requested from: Andrea Piccolo, Boris Kaus
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
- Introduced by (git blame at HEAD):
  - Garrett Ito — 11be8a9e Locate_dike: nPtr set in routine... (2023-04-30) — x loop bound uses sx+ny instead of sx+nx
  - Garrett Ito — e58ddd84 dike.cpp 2D elliptical Gaussian weighted filter (2023-04-30) — wrote time-averaging block
- At HEAD: still present, src/dike.cpp:1331. Note: Same author both commits.
- Feedback requested from: Garrett Ito
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
- Introduced by (git blame at HEAD):
  - Garrett Ito — f1da6107 Wrote code for Behn & Ito Dike heat terms (2021-07-27); d9d3a717 Compiling somehow stuck on parsing.h (2021-07-27) — introduced k = kfac*k with kfac summing only dike phRat
  - Jana Schierjott — b7c8c459 dike heating dependent on varying M (2021-11-18); 7fe091f2 removed print statements (2021-09-01) — reworked kfac accumulation inside phase loop; k=kfac*k per dike
- At HEAD: still present, src/dike.cpp:463-479. Note: Behn & Ito (2008) formulation originally Garrett's; Jana extended it to multiple dikes/M.
- Feedback requested from: Garrett Ito
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
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 348d9a2e passive tracer changes, linear interpolation (2020-11-27) — wrote pressure interpolation with vz pattern (K, ncz)
  - Andrea Piccolo — beeec8c5 melt fraction computation on passive tracers (2020-11-10) — wrote temperature interpolation with vz pattern
  - Anton Popov — 076ba27c Get rid of all global buffers (2026-07-01) — renamed variables in same lines; bug preserved
- At HEAD: still present, src/passive_tracer.cpp:560-561. Note: Anton's 2026 refactor is blame target but logic is Andrea's.
- Feedback requested from: Anton Popov, Andrea Piccolo
- Status: fixed (group A)

### S-025 [major, confirmed] Observation points in the same cell overwrite each other's weights
- Tags: major, confirmed, logic, quick
- Where: src/adjoint.cpp:2342-2349, 2387-2394, 2433-2440 (and `xini` lines below each)
- What: projection weights are set with `=` in the local array. Two points that share nodes
  (same or neighbouring cell on one CPU) overwrite each other; only the last point counts there.
- Why it matters: with many points (e.g. dense surface data) part of the data is silently lost
  from the cost and the gradient.
- Fix idea: `+=` (the local-to-global step already uses ADD_VALUES), or one projection per point.
- Introduced by (git blame at HEAD):
  - Georg Reuber — 2da29341 Initial creation of the branch with first working version (2020-05-08); cf0ea953 Added the adjoint gradient computation and the TAO framework (2018-01-11) — Wrote the point-projection weight assignment with '=' in AdjointPointInPro (rewritten 2020, same logic).
- At HEAD: still present, src/adjoint.cpp:2348-2355, 2393-2400, 2439-2446 (llproX/Y/Z '=' assignments; xini copies 2370-2377, 2415-2422, 2461-2468). Note: 2da29341 is a substantive rewrite by the same author, not a bulk reformat.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: please throw an error if this happens; don't think this was ever an issue. Summing may not be trivial withn adjoint.

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
- Introduced by (git blame at HEAD):
  - Boris Kaus — 05bc932d Clarified initial overview of parameters; added routines to add/delete parameters from option database (2020-05-14); 3a9e100d WIP: TAO (1); ensured that parameters are set on the command-line for TAO (2020-05-17) — Introduced 'if(ts)' bound test with InitGuess defaults; 3a9e100d added Ab=1 trigger for BLMVM.
- At HEAD: still present, src/adjoint.cpp:633-661 (if(ts) bounds, else Ubar/Lbar=par_val; Ab=1 at 640, 655); used at 1049-1058. Note: Before 05bc932d Ubar/Lbar were stored unconditionally (default 0); the default-to-InitGuess logic is new there.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: try the fix suggested (not sure if it'll converge with INF)

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
- Introduced by (git blame at HEAD):
  - Georg Reuber — 2da29341 Initial creation of the branch with first working version (2020-05-08); cf0ea953 Added the adjoint gradient computation and the TAO framework (2018-01-11); fdbed43e Clean up of LaMEM.cpp. (2018-01-15) — Changed fcconv size from 1501 to maxit (default 50); original writes and 1500 print-only check.
  - Boris Kaus — f38129dc fixed a few memory leaks (2020-05-18) — Moved fcconv creation into AdjointVectorsCreate; size unchanged.
- At HEAD: still present, src/adjoint.cpp:1159 (fcconv size maxit), 1438 (descent write), 1534 (TAO write), 1456/1539 (count checks), 1094-1099 (print). Note: Out-of-bounds arose when 2da29341 replaced the fixed 1501 entries by maxit; TAO still counts to 1500.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: fix

### S-028 [major, confirmed] `factor2array[51]` is indexed up to 100 parameters
- Tags: major, confirmed, memory, quick
- Where: src/objFunct.h:85; src/adjoint.cpp:1294, 1325, 1346, 1408, 1449
- What: the line-search factors are a fixed array of 51, but loops run to `mdN` (up to
  `_MAX_PAR_ = 100`).
- Why it matters: with more than 51 parameters and `Inversion_EmployTAO = 0`, writes past the
  array into the next `ModParam` fields (`maxfac`, `vel_scale`, ...).
- Fix idea: size it `_MAX_PAR_`.
- Introduced by (git blame at HEAD):
  - Georg Reuber — 2da29341 Initial creation of the branch with first working version (2020-05-08); 3492efcd fixing a bunch of things in the GD algorithm (2020-06-08) — Introduced fixed factor2array[51] and per-parameter loops; 3492efcd rewrote the loops using it.
- At HEAD: still present, src/objFunct.h:85 (factor2array[51]); uses at src/adjoint.cpp:1296, 1328, 1343, 1356, 1415, 1445. Note: Blame on objFunct.h:85 shows Boris Kaus b49a6e9f (2020-05-22, WIP scaling) but that only reflowed comments.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: implement suggested fix

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
- Introduced by (git blame at HEAD):
  - Anton Popov — 94177af5 . (2017-02-16); 4e37814a Added output of overpresure field (2016-08-26); d9bebfb1 ... (2017-01-13) — Wrote lp_lith/lp_pore computation (owned cells only, zero ghosts) and pore-pressure edge averaging.
  - Jianfeng Yang ("Jeff") — 5157b028 Compute lithostatic pressure; pressure bound for Drucker-Prager plasticity (2016-10-11) — Introduced 4-cell edge averaging of lithostatic pressure (incl. outside ghosts) for plasticity limit.
- At HEAD: still present, src/JacResAux.cpp:244-321 (lith), 339-411 (pore); edge averages src/JacRes.cpp:1195-1198, 1302-1305, 1410-1413. Note: Jeff's commit (JacRes.c then) created the edge averaging; Boris Kaus bd65f30e (2016) added rho_lithos viscosity use but not the averaging.
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — af6dfe3d Anton allowed me to push the adiabatic heating into the main (2019-07-12); acc54759 I apologize. I forgot to plug Ha in the thermal equation (2019-07-18) — Added adiabatic-heating pressure-gradient term using neighbour-cell lithostatic pressure.
- At HEAD: still present, src/JacResTemp.cpp:517-526 (bdpdx.. lines 520-522 use P[...] incl. ghost neighbours). Note: Depends on zero outside ghosts of lp_lith (see S-029, Anton Popov / Jeff).
- Feedback requested from: Andrea Piccolo
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
- Introduced by (git blame at HEAD):
  - Georg Reuber — 61aa7159 Added phase diagrams to LaMEM (2018-01-19); e9b938f9 Small bugfix in phase diagrams (2018-01-23) — Wrote the (indP-1) row indexing, index<1 and upper clamps in the interpolation.
  - Boris Kaus — ddb7f6aa rewriting part of the routines for code clarity (2018-01-19) — Rewrote variable access (nT/nP, minP/dP); kept Georg's index logic.
- At HEAD: still present, src/constEq.cpp:1008-1060 (setDataPhaseDiagram; ind[] at 1047-1050, clamps 1019-1058). Note: ddb7f6aa blame is a same-day rewrite on top of 61aa7159; the (indP-1) logic is in 61aa7159 itself.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: that is imndeed an issue; implement suggested fix. Will change the test results

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
- Introduced by (git blame at HEAD):
  - Wenrong Cao — 749ae824 Add slope-dependent erosion to the free surface (v3.1.0) (2026-08-13) — Implemented FreeSurfAppSlopeErosion with centered-difference slope.
- At HEAD: still present, src/surf.cpp:1229-1253 (centered gradient + erosion); slope scan 1173-1188. Note: New feature (single commit); function skeleton shares loops with Wenrong Cao 1a44f790 topographic diffusion.
- Feedback requested from: Wenrong Cao
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
- Introduced by (git blame at HEAD):
  - Wenrong Cao — 749ae824 Add slope-dependent erosion to the free surface (v3.1.0) (2026-08-13) — Wrote the sub-step rule limiting vertical drop to half the grid spacing.
- At HEAD: still present, src/surf.cpp:1199-1211 (E_max, dt_cfl = 0.5*d_min/E_max, nsolve)
- Feedback requested from: Wenrong Cao
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
- Introduced by (git blame at HEAD):
  - Yun Luo — e4778c61 New coupling module based on LaMEM 3.0.1 (2026-08-02) — Added FastScape coupling: smoothing after FastScapeRun while FastScape restarts from own saved topography.
- At HEAD: still present, src/LaMEMLib.cpp:749-752 (FastScapeRun then FreeSurfSmoothMaxAngle); src/fastscape.cpp:1590-1606 (gtopo_fs reused). Note: FreeSurfSmoothMaxAngle itself is old (Anton Popov 2016); the conflict arises only from the new coupling.
- Feedback requested from: Yun Luo
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
- Introduced by (git blame at HEAD):
  - Anton Popov — 8763ccaf Multiple updates from concurrent branches (2017-01-11) — Split delete path into AVDDeletePointsMV, copying area=tclaimed without the claim-counting loop.
  - Adina Püsök — 9edd2190 bug fix for Random noise that was making the code hanging (2015-04-01); 766f54c3 1. Marker control (solves the Empty Control Volume error) (2015-03-30) — Wrote original area/tclaimed sort in combined AVDInjectDeletePoints (where counting was present).
- At HEAD: still present, src/AVD.cpp:1186-1214 (AVDDeletePointsMV; area[i]=tclaimed at 1195). Note: git blame -C attributes the line to Adina Pusok (copied); bug (tclaimed never counted in delete path) arose in Anton's 2017 split.
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Arne Spang — d99dc47b Removed a typo in a comment and made scaling operation (2025-02-24); 9caf07c6 Updated confusing naming of phase transitions (2025-02-24) — Changed slope scaling from 1e6*(T/stress_si) to /(stress/temperature), wrong in SI units.
  - Andrea Piccolo — 1bc3ab4d Testing the new boundary conditions, still working on it (2020-10-04); a97d96fc Return to a marker based approach (2020-07-03) — Added T0 + Tshift conversion and predefined laws with T0 in deg C.
- At HEAD: still present, src/phase_transition.cpp:592-594 (scaling), 667-711 (predefined laws in deg C / MPa/K); scaling.cpp. Note: Slope bug is Arne's 2025 edit (earlier formula was SI-aware); Tshift/T0-in-degC issue dates to Piccolo 2020.
- Feedback requested from: Andrea Piccolo, Anton Popov
- Status: open

### S-037 [major, confirmed] NotInAirBox uses only the first segment's `zbounds`
- Tags: major, confirmed, logic, medium
- Where: src/phase_transition.cpp:1303, 1322-1323, 1334 (commit 6124aa7d)
- What: `zbounds` is read for every segment (`2*nsegs` values) but the inside test and the
  temperature profiles always use `zbounds[0]` and `zbounds[1]`.
- Why it matters: segments with a different depth range (as in the documented example,
  input_file.dat:781) get the wrong vertical extent.
- Fix idea: store per-cell z bounds like `celly_xboundL/R`, or document that only one z range
  is supported and read only two values.
- Introduced by (git blame at HEAD):
  - Garrett Ito — c058f984 phase_transition: working with old input file using ph->bounds (2022-06-25); dea561da phase_transition.cpp; ph->celly_xboundL and xboundR kept (2022-06-26) — Made xbounds per segment and zbounds read 2*nsegs but kept using zbounds[0]/[1].
  - Boris Kaus — d7142fe3 Finished implementation of Box-like phase transition (2020-11-17) — Original single-box NotInAirBox with temperature profiles using one z range.
- At HEAD: still present, src/phase_transition.cpp:1302-1304 (zbounds[0],[1] test), 1322-1323, 1334; read at 413. Note: Bug appeared when Garrett added multi-segment support but not per-segment z bounds.
- Feedback requested from: Garrett Ito, Boris Kaus
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
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — dd3b7cf1 I commit the first stage of the work. Still need to t... (2020-10-21); 7968d4f4 Update with the old files and start working on it (2020-10-21) — Placed ADVAdvectPassiveTracer before ADVRemap and used markstart/markind in tracer code.
- At HEAD: still present, src/LaMEMLib.cpp:720-721 (call before ADVRemap); src/passive_tracer.cpp:585-607, 1032-1049, 904-938. Note: Call order unchanged since the introduction of passive tracers.
- Feedback requested from: Andrea Piccolo
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
- Introduced by (git blame at HEAD):
  - Anton Popov — b7f5b864 default solver settings (update) (2025-08-22); 66d86f44 ... (2025-08-25); 87f79983 Default solver settings complete! (2025-08-26) — Defaults now create Scaling/FDSTAG and read act_temp_diff from input file, ignoring command line.
  - Boris Kaus — ddb305fd added default solver options, which can now be specified (2018-01-16) — Original clear-options / set-defaults / reinsert-command-line structure.
- At HEAD: still present, src/options.cpp:937-960 (setSolverOptions clear -> defaults -> reinsert); defaults building scaling/grid at 120-165, 251-256. Note: Structure is Boris's 2018 design; the grid-dependence that makes it matter is Anton's 2025 code.
- Feedback requested from: Anton Popov, Boris Kaus
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
- Introduced by (git blame at HEAD):
  - Georg Reuber — 2da29341 Initial creation of the branch with first working version (2020-05-08); cf0ea953 Added the adjoint gradient computation and the TAO framework (2018-01-11) — Wrote the scaled-xini point-velocity cost and dF/dx = w*(x-xini).
- At HEAD: still present, src/adjoint.cpp:2357-2380 (xini = Ae/f, f=v/v_node; Vy/Vz at 2402-2426, 2448-2472), 1650-1676
- Feedback requested from: Boris Kaus
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
- Introduced by (git blame at HEAD):
  - Boris Kaus — 3e17641e Generalized the AdjointGet_F_dFdu_Center routine (2020-06-05); 9f902278 Bugfix: Gradients for PSD are now correctly computed (2020-06-05); 0215c13f added scaling for Adjoint/Strainrates (2020-06-08) — Strain-rate target stored unscaled; incomplete (non-deviatoric) derivatives; scaling added only to printed value.
  - Georg Reuber — 625d00d2 added FD kernel calculation for eta and n (2020-06-06); 67738f08 Added PSD angle derivative calculation (2020-06-02) — Added the strain-rate observation framework and E2nd branch.
- At HEAD: still present, src/adjoint.cpp:868-875 (Ae = ts unscaled), 3455-3582 (derivatives, E2nd 'not yet implemented'); src/JacRes.cpp:675-677. Note: Effective strain-rate with elastic term (JacResGetEffStrainRate) is Anton Popov's JacRes design.
- Feedback requested from: Boris Kaus
- Status: open
- Boeis: I a, niot sure we ever used strainrate for inversion, so could be that this is wrong. If velocity is propoerly non-dimensionalized., this should be fixed 

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
- Introduced by (git blame at HEAD):
  - Georg Reuber — cf0ea953 Added the adjoint gradient computation and the TAO framework (2018-01-11) — Originally wrote Ax/Ay/Az += v*dt inside the on-processor block.
  - Boris Kaus — d2cb841f changes to make a GenericInversion work (2020-06-19) — Restructured block (current blame owner); behaviour unchanged.
- At HEAD: still present, src/adjoint.cpp:1744-1749 (inside lrank == 13 block). Note: Blame shows d2cb841f but -S shows the same expression in cf0ea953.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: I don't think we have tested adjoint in parallel; yes this cpu;ld nbe wrom. Perhaps fix it (low priority)

### S-043 [major, possible] FastScape grid doesn't follow a stretching LaMEM box
- Tags: major, possible, physics, design
- Where: src/fastscape.cpp:458-469 (ranges set once), 646 (mode 0 "update range" never called)
- What: rangeX/Y/Z and the FastScape grid are set once at start-up. With background strain
  (BCStretchGrid), LaMEM's box changes, but topography is still mapped node by node and
  clamped to the old z range. The velocities sent to FastScape also include the stretching part.
- Why it matters: with bg strain + FastScape, topography ends up at the wrong x/y position and
  is clamped to wrong limits. Not tested; t37 has no background strain.
- Fix idea: update the range each step (mode 0), or reject bg strain with surf_mode = 2.
- Introduced by (git blame at HEAD):
  - Yun Luo — e4778c61 New coupling module based on LaMEM 3.0.1 (2026-08-02) — Wrote the FastScape grid setup with fixed ranges and unused update-range mode.
- At HEAD: still present, src/fastscape.cpp:456-469 (ranges set once), 646-668 (CreateSurfaceGrid mode 0 never called; only modes 1 and 2 at 178, 2948)
- Feedback requested from: Yun Luo
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
- Introduced by (git blame at HEAD):
  - Anton Popov — 1f0d3a8b update (2017-05-19); 227fbb31 update (2017-05-12); 139c4a1a constrained cells form parallel files (2017-08-02) — Wrote the fixed-phase/fixed-cell loops setting i+1/j+1/k+1 faces.
- At HEAD: still present, src/bc.cpp:2594-2606 (BCApplyPhase), 2647-2659 (BCApplyCells)
- Feedback requested from: Anton Popov
- Status: open

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
- Introduced by (git blame at HEAD):
  - Anton Popov — 80a54602 PBC update (2025-09-03); 91006e77 PBC update (2025-08-30); eb35f055 Split periodic advection and periodic grid flags (2025-08-30) — Set actx->periodic for Exy and forbade periodic grid together with Exy.
- At HEAD: still present, src/advect.cpp:277-280 (periodic flag with ExyNumPeriods), 1084-1090; src/bc.cpp:904-909. Note: Arose in the 2025 PBC refactor which replaced the earlier cycle_adv mechanism (7b0ccd0e, 2022, Anton Popov).
- Feedback requested from: Anton Popov
- Status: open

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
- Introduced by (git blame at HEAD):
  - Anton Popov — 9a108627 Explicit time stepping + elasticity + kinematic bc block (2017-03-14) — Introduced immediate dt update to dt_next before marker advection when dt is not fixed.
- At HEAD: still present, src/tssolve.cpp:307-308 (if(!fix_dt) ts->dt = dt_next); call at src/LaMEMLib.cpp:697
- Feedback requested from: Anton Popov
- Status: open

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
- Introduced by (git blame at HEAD):
  - Georg Reuber — cf0ea953 Added the adjoint gradient computation and the TAO framework (2018-01-11); 2da29341 Initial creation of the branch with first working version (2020-05-08) — Wrote surface-velocity observation loops over corner indices and copied FindPointInCell with SETERRQ.
- At HEAD: still present, src/adjoint.cpp:2505-2676 (Ap == 3 block; K lookup 2548, 2605, 2661); FindPointInCellAdjoint 102-125
- Feedback requested from: Boris Kaus
- Status: open
- Boris: would be nice to fix

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
- Introduced by (git blame at HEAD):
  - Georg Reuber — f34e7d65 Clean up of creation/destruction (1) (2020-05-18); 2da29341 Initial creation of the branch with first working version (2020-05-08); fdbed43e Clean up of LaMEM.cpp. (2018-01-15) — Created parameter vectors with fixed per-CPU local size _MAX_PAR_ (and mdI for obs vectors).
- At HEAD: still present, src/adjoint.cpp:1152-1156 (Lb, Ub, val, P, grad with _MAX_PAR_ local size); 1186-1189 (aop->vx..sty with mdI local size). Note: Boris Kaus f38129dc (2020-05-18) only relocated creation.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: would be nice to fix

## Minor

### M-001 [minor, confirmed] Adiabatic heating efficiency is applied twice
- Tags: minor, confirmed, logic, quick
- Where: src/JacResTemp.cpp:525 and :534
- What: `Ha = AdiabHeat*(T*alpha*...)`, then `Ha = AdiabHeat*Ha`. The output value
  (`svBulk->Ha`) has the factor once, but the energy equation gets it squared.
- Why it matters: no effect for 0 or 1. With `Adiabatic_Heat = 0.5`, the model uses 25%.
- Fix idea: remove line 534.
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 268eddb0 A little bug that need to be checked. (2021-01-27) — added second Ha = AdiabHeat*Ha (line 534) after svBulk->Ha stored
  - Andrea Piccolo — acc54759 I forgot to plug Ha in the thermal equation... (2019-07-18) — first applied AdiabHeat factor inside Ha expression (line 525)
- At HEAD: still present, src/JacResTemp.cpp:525,534. Note: Both commits by Piccolo; the double application arose in 268eddb0 when line 534 was added.
- Feedback requested from: Andrea Piccolo
- Status: open

### M-002 [minor, confirmed] Top pressure shift isn't area-weighted and overwrites `p_shift`
- Tags: minor, confirmed, physics, quick
- Where: src/JacRes.cpp:558-602
- What: the top cell layer is averaged with equal weights. With `act_p_shift` on, the result
  replaces any user-given `p_shift` without a warning.
- Why it matters: on grids with variable spacing, the mean top pressure isn't really zero.
  The user's input is silently ignored.
- Introduced by (git blame at HEAD):
  - Boris Kaus — ed8eab85 reactivated the pressure shift... (2020-11-08) — equal-weight top-layer pressure average; overwrites ctrl.pShift
- At HEAD: still present, src/JacRes.cpp:565-600 (JacResGetPShift, lines 582,599). Note: Later commits (Popov 076ba27c) only reformatted/PetscCall.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: indeed a bug; please fix

### M-003 [minor, confirmed] `Phase_Melt` density ignores the `mfmax` cap
- Tags: minor, confirmed, logic, quick
- Where: src/constEq.cpp:734-742
- What: in the `Phase_Diagram_melt` branch, `mf` is capped by `mfmax` but then not used.
  The density uses `Pd->mf` directly.
- Fix idea: use `mf` on line 741.
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — da9950e0 Cleaning a bit the branch, and small bug fix in the consteq.cpp (2021-01-13) — added mf cap copy but kept Pd->mf in melt-density line 741
- At HEAD: still present, src/constEq.cpp:739-741 (Phase_Diagram_melt branch). Note: Cap pattern copied from Kaus adde8dc1 (2020-11-10) for the non-melt branch; melt branch line never switched to mf.
- Feedback requested from: Andrea Piccolo
- Status: open

### M-004 [minor, confirmed] Healing defaults set only when printing is on
- Tags: minor, confirmed, logic, quick
- Where: src/phase.cpp:223-232
- What: `APSheal2 = APS2` and `healTau2 = healTau` are set inside `if(PrintOutput)`. The
  adjoint re-reads the database with `PETSC_FALSE` (adjoint.cpp:2881), so that copy has no
  second healing stage.
- Fix idea: move the default assignment out of the print block.
- Introduced by (git blame at HEAD):
  - Garrett Ito — c2f25559 Updates since last sync with bitbucket/lamem kaus ... (2023-08-06) — set APSheal2/healTau2 defaults only inside the PrintOutput branch
- At HEAD: still present, src/phase.cpp:215-220 (inside if(PrintOutput)). Note: Adjoint re-read adjoint.cpp:2881 uses PETSC_FALSE so defaults are skipped; confirmed still present.
- Feedback requested from: Garrett Ito
- Status: open

### M-005 [minor, confirmed] Sediment model 3 assumes the box starts at x = 0
- Tags: minor, confirmed, logic, quick
- Where: src/surf.cpp:1531
- What: `dz_x = -dz/BoxWidth*x + dz1` uses x, not x - bx. The rate is rate1 at the left edge
  and rate2 at the right edge only if bx = 0. The comment above it (Gemmer Gaussian) is
  copied from model 2 and is wrong: this is a linear ramp. Unlike model 1, there's no level check.
- Fix idea: use (x - bx).
- Introduced by (git blame at HEAD):
  - Tobias Baumann — e380d6ae typo (2019-04-04) — wrote dz_x = -dz/BoxWidth*x + dz1 using x, not x-bx
  - Tobias Baumann — e87ce628 added another sedimentation mechanism (differential loading) (2019-04-04) — introduced sediment model 3 and copied Gemmer comment
- At HEAD: still present, src/surf.cpp:1531 (SedimentModel==3). Note: tosaba/baumann@uni-mainz.de mapped to Tobias Baumann (name inferred, not seen in git).
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-006 [minor, confirmed] Initial surface perturbation options have unit and NaN problems
- Tags: minor, confirmed, input, parallel, quick
- Where: src/surf.cpp:1711-1754
- What: `-FreeSurf_Wavelength` and `-FreeSurf_AmplNoise` are not scaled (only AmplCos is
  divided by the length scale), so they're wrong unless the length scale is 1. With noise only
  (wavelength 0), 2*pi/0 gives inf, cos(inf) = NaN, and 0*NaN = NaN in every node. Each rank uses
  the same seed, so the noise pattern repeats on every CPU tile and depends on CPU count.
- Fix idea: scale all three, skip the cosine if wavelength is 0, and seed by global node index.
- Introduced by (git blame at HEAD):
  - Boris Kaus — 01b07ee4 added a Rayleigh-Taylor Instability benchmark to the testing suite (2020-05-26) — added -FreeSurf_Wavelength/AmplCos/AmplNoise options; only AmplCos scaled
  - Boris Kaus — ada81a7a added a reference figure for RTI test; same-seed random noise (2020-05-26) — added identical-seed noise on every rank
- At HEAD: still present, src/surf.cpp:1711-1754 (FreeSurfSetInitialPerturbation). Note: Popov 745ab1bb (2025-11-02) only touched the line for periodic BC.
- Feedback requested from: Boris Kaus
- Status: open
- BNoris: please fix

### M-007 [minor, confirmed] Topographic diffusion boundary is not zero-flux
- Tags: minor, confirmed, physics, quick
- Where: src/surf.cpp:1566-1567, 1651, 1662
- What: boundary nodes get no x (or y) term, but their neighbours still exchange with them. That
  isn't a Neumann condition (which would need a half-cell term at the boundary node), and the
  total volume isn't kept. The comment says "zero-flux Neumann".
- Fix idea: use a ghost node with z[-1] = z[1], or fix the comment.
- Introduced by (git blame at HEAD):
  - Wenrong Cao — 1a44f790 Add topographic diffusion feature to free surface module (2026-04-04) — boundary nodes skip x/y flux terms; comment claims zero-flux Neumann
- At HEAD: still present, src/surf.cpp:1566-1567,1649-1669 (FreeSurfAppTopoDiffusion)
- Feedback requested from: Wenrong Cao
- Status: open

### M-008 [minor, confirmed] Erosion and sedimentation levels overshoot by up to one step
- Tags: minor, confirmed, logic, quick
- Where: src/surf.cpp:1004-1008, 1068-1072, 1346-1350
- What: a node just above the erosion level is lowered by the full rate*dt, so it can end up to
  rate*dt below the level. Sedimentation does the same thing in the other direction.
- Fix idea: clamp to the level, e.g. z = max(z - dz, level).
- Introduced by (git blame at HEAD):
  - Xinxin Wang — c0faa50d added the described erosion/sedimentation levels to the surface processes (2021-03-24) — added 'if(z > level) z -= dz' / 'if(z < level) z += dz' without clamping to level
  - Wenrong Cao — 843e9eb2 Add spatially-limited erosion model (erosion_model = 3) (2026-01-26) — copied same unclamped level test into erosion model 3
  - Tobias Baumann — c645c37b so far we only had infinitely fast erosion; added simple erosion at constant rate (2019-06-18) — original uniform rate*dt advection z -= dz
- At HEAD: still present, src/surf.cpp:1004-1008, 1068-1072, 1346-1350
- Feedback requested from: Wenrong Cao, Anton Popov
- Status: fixed (group A)

### M-009 [minor, possible] Surface outside the box silently gets zero velocity
- Tags: minor, possible, logic, quick
- Where: src/surf.cpp:461 (`z >= bz && z < ez`)
- What: a surface node at or above the top of the box (or below the bottom) isn't found on any
  rank, so its velocity is 0. FreeSurfAdvectTopo doesn't clamp, and FastScape clamps to exactly
  the top (fastscape.cpp:1957, 2093), which is excluded by `z < ez`.
- Fix idea: clamp z into [bz, ez] on the top/bottom ranks before the lookup, or warn.
- Introduced by (git blame at HEAD):
  - Anton Popov — 157753dd Bug fixes. (2017-03-24) — introduced half-open test z >= bz && z < ez
- At HEAD: still present, src/surf.cpp:461 (FreeSurfAdvectTopo velocity lookup). Note: FastScape clamp to exactly top (Yun Luo module) feeds this edge case.
- Feedback requested from: Anton Popov
- Status: fixed (group C, low value)

### M-010 [minor, confirmed] Memory leaks in FastScape code
- Tags: minor, confirmed, memory, quick
- Where: src/fastscape.cpp:1210, 1234-1242 vs 1310-1317; 3100-3110
- What: VecScatterCreateToZero makes the output Vecs on every rank, but they're destroyed only on
  rank 0, so each step leaks 3 Vecs per non-root rank (4 at step 0). FastScapeDestroy frees only
  one of `ncoor_ori_x` / `ncoor_ori_y`, but 3D non-uniform grids allocate both.
- Introduced by (git blame at HEAD):
  - Yun Luo — e4778c61 New coupling module based on LaMEM 3.0.1 (2026-08-02) — added FastScape coupling file; Vecs from VecScatterCreateToZero leaked off-root; one ncoor_ori freed
- At HEAD: still present, src/fastscape.cpp:1210,1234-1242,1305-1317; 3100-3110 (FastScapeDestroy). Note: Whole file arrived in e4778c61; later Kaus/Popov commits are formatting/restructuring.
- Feedback requested from: Yun Luo
- Status: open

### M-011 [minor, possible] Non-uniform FastScape interpolation can skip a node
- Tags: minor, possible, logic, quick
- Where: src/fastscape.cpp:1055-1079
- What: if FindIndexForInterpolation misses a point (e.g. the last node, a bit outside the range
  due to roundoff), that node isn't written, so it keeps its old value (0 at the start).
- Fix idea: clamp coordinates into the range, or report an error.
- Introduced by (git blame at HEAD):
  - Yun Luo — e4778c61 New coupling module based on LaMEM 3.0.1 (2026-08-02) — interpolation loop skips node when FindIndexForInterpolation misses
- At HEAD: still present, src/fastscape.cpp:1046-1081 (if(indices.found) without else)
- Feedback requested from: Yun Luo
- Status: open

### M-012 [minor, confirmed] `interpStretch` can leave stretch factors unset after the last control polygon
- Tags: minor, confirmed, logic, memory, quick
- Where: src/tools.cpp:628 (used by marker.cpp:1696)
- What: the check is `CtrlPoly[last] != numCtrlPoly-1` but should be `!= numPoly-1`. If the last
  control polygon happens to sit at index numCtrlPoly-1 (e.g. control polygons 1 and 2 of a
  10-slice volume), slices after it get no Sx/Sy. Those are uninitialized stack values that
  stretchPolygon then uses.
- Fix idea: compare with numPoly-1, or simply always fill the tail.
- Introduced by (git blame at HEAD):
  - Arne Spang — a88df698 added geometry manipulation (2020-07-08) — wrote CtrlPoly[last] != numCtrlPoly-1 (should be numPoly-1)
- At HEAD: still present, src/tools.cpp:628 (interpStretch)
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-013 [minor, confirmed] Phase 0 can't be used as permeable / plume mantle inflow phase
- Tags: minor, confirmed, input, quick
- Where: src/bc.cpp:848
- What: the check `!bc->phase_inflow_bot` treats phase 0 as "not set", because the default
  is 0 (zeroed struct).
- Why it matters: models where the mantle is phase 0 stop with a wrong error.
- Fix idea: default to -1 and test `< 0`.
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 6d86dac0 small correction (2021-02-09) — check '!bc->phase_inflow_bot' treats phase 0 as unset
  - Andrea Piccolo — da9950e0 Cleaning a bit the branch... (2021-01-13) — added permeable_phase_inflow option with default -1
- At HEAD: still present, src/bc.cpp:848
- Feedback requested from: Andrea Piccolo
- Status: open

### M-014 [minor, confirmed] `relax_dist = 0` with `bvel_face_out` gives NaN at top/bot
- Tags: minor, confirmed, logic, quick
- Where: src/bc.cpp:2097-2099, 2161-2163
- What: with the default `relax_dist = 0`, a cell center with `z == top` (or `z == bot`)
  passes the ramp test and computes `velin - (velin/0)*0 = NaN`.
- Why it matters: rare, but grids with round coordinates can hit it exactly.
- Fix idea: skip the ramp when `relax_dist <= 0`.
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 1bc3ab4d Testing the new boundary conditions, still working on it (2020-10-04) — wrote velin-(velin/relax_dist)*(z-top) ramp without relax_dist>0 guard
- At HEAD: still present, src/bc.cpp:2097-2098, 2161-2162
- Feedback requested from: Andrea Piccolo
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
- Introduced by (git blame at HEAD):
  - Panupong Kongpet — 50a251b3 Add 3D Guassian temperature in bc (2026-01-22) — 3D BC uses exp(-r2/(2R2)), differs from marker exp(-r2/R2); Tbot used unset
  - Andrea Piccolo — 91955c88 Other small correction on the marker injection routine (2021-01-22) — marker T_inflow = Tbot+(Tplume-Tbot)*exp(-r2/R2) with Tbot possibly -1
- At HEAD: still present, src/bc.cpp:1507,1519 (BC) vs 3005,3017 (marker injection). Note: Ito c39179e2 (2026-02-12) later changed only the 3D radius cutoff.
- Feedback requested from: Andrea Piccolo
- Status: open

### M-016 [minor, confirmed] Typos in BC options are silently ignored; `exy` drops no-slip on vx
- Tags: minor, confirmed, input, quick
- Where: src/bc.cpp:675-679, 712-713, 1784-1802
- What: an unknown `bvel_face` (e.g. "left") leaves `face = 0`, so no inflow BC is applied
  and no message. Same for `bvel_temperature_inflow`. With `exy_strain_rates` set, the
  `else` branch that applies no-slip on vx at front/back/bottom/top is skipped, so those
  no-slip flags are silently ignored for vx.
- Fix idea: error on unknown strings; warn or error when `noslip` and `exy` are combined.
- Introduced by (git blame at HEAD):
  - Anton Popov — 94ea1409 Tangential boundary velocity support (2024-05-07) — added if(Exy) ... else so no-slip vx BCs skipped when exy set
  - Boris Kaus — c1ecd9d3 changed printing statement with inflow BC's... (2020-11-17) — bvel_face strcmp chain; unknown string leaves face=0
  - Andrea Piccolo — f460e203 Fixed thermal inflow bug, change input file (2020-11-13) — bvel_temperature_inflow strcmp chain; typo silently ignored
- At HEAD: still present, src/bc.cpp:675-679, 712-713, 1775-1802
- Feedback requested from: Anton Popov, Andrea Piccolo, Boris Kaus
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
- Introduced by (git blame at HEAD):
  - Anton Popov — 38439615 Viscosity coarsening (done) (2025-04-19) — wrote 'if(I > mnx) I--' edge clamp (should be >=) in rediscretized coarsening
- At HEAD: still present, src/matData.cpp:497-498, 529-530, 561-562 (MatDataRestrictParam3D)
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-018 [minor, confirmed] Some SNES/KSP failure reasons are not printed or not treated as severe
- Tags: minor, confirmed, logic, quick
- Where: src/nlsolve.cpp:430-515 (SNESPrintConvergedReason)
- What: `SNES_DIVERGED_DTOL`, `SNES_DIVERGED_TR_DELTA` and `SNES_DIVERGED_JACOBIAN_DOMAIN`
  print the "FAILED" header but no reason. Under `SNES_DIVERGED_LINEAR_SOLVE`,
  `KSP_DIVERGED_PC_FAILED` (e.g. zero pivot in LU, failed coarse solve) is not in the severe
  list, so with `-snes_continue_on_fail` the run goes on with a broken preconditioner.
- Fix idea: add the missing cases; treat `KSP_DIVERGED_PC_FAILED` as severe.
- Introduced by (git blame at HEAD):
  - Anton Popov — 98e959fc stop in case that severe SNES and/or KSP divergence is detected (2019-03-29) — severe-KSP list lacks KSP_DIVERGED_PC_FAILED
  - Anton Popov — 88b60a16 code integration & functionality restore update version 0.0 (2014-09-09) — original reason chain without DTOL/TR_DELTA/JACOBIAN_DOMAIN cases
- At HEAD: still present, src/nlsolve.cpp:430-515 (SNESPrintConvergedReason)
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Anton Popov — 324fb7bf update (2024-10-11) — MFFD operator calls SNESComputeFunction (counted by SNES)
  - Anton Popov — b7f5b864 default solver settings (update) (2025-08-22) — option set '-snes_max_funcs unlimited' (mitigates)
- At HEAD: partly mitigated, src/nlsolve.cpp:284 (MatMFFDSetFunction) ; mitigation src/options.cpp:59. Note: MFFD still counts calls, but options.cpp now sets -snes_max_funcs unlimited, so the cap is not hit by default.
- Feedback requested from: Anton Popov
- Status: won't fix (group C): mitigated by `-snes_max_funcs unlimited`

### M-020 [minor, confirmed] ADVCheckCorners reads an empty buffer
- Tags: minor, confirmed, memory, quick
- Where: advect.cpp:1485, 1526-1539 (commit 6124aa7d)
- What: if the cell and all its local neighbours have no markers, `n == 0`, the parent loop
  never runs, and `markers[sind]` is read from a zero-size buffer with an unset `sind`.
- Why it matters: crash or garbage marker in emptied in/outflow regions, more likely at rank
  corners.
- Fix idea: skip the corner (or use a default phase) when n == 0.
- Introduced by (git blame at HEAD):
  - Adina Püsök — e5016a08 added local neighbor cells influence for corner marker injection (2015-06-15) — neighbour-cell loading; n can be 0 so markers[sind] reads empty buffer
  - Adina Püsök — 766f54c3 Marker control (solves the Empty Control Volume error) (2015-03-30) — original corner injection with markers[sind], sind unset when no markers
- At HEAD: still present, src/advect.cpp:1539 (markers[sind], ADVCheckCorners)
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-021 [minor, confirmed] AVD never refills empty volumes
- Tags: minor, confirmed, logic, quick
- Where: AVD.cpp:753-754, 1136-1137 (commit 6124aa7d)
- What: `ninj += min(nmin-n, n)`: with n = 0 nothing is injected, and at most n per pass.
- Why it matters: under `mctrl = avd` an empty volume stays empty, with no warning.
- Introduced by (git blame at HEAD):
  - Anton Popov — 8763ccaf Multiple updates from concurrent branches (2017-01-11) — wrote 'if((nmin-n)>n) ninj += n' logic injecting at most n markers
- At HEAD: still present, src/AVD.cpp:753-754. Note: Audit's second location (AVD.cpp:1136) no longer has ninj+= at HEAD.
- Feedback requested from: Anton Popov
- Status: not an issue (AVD batch): no parent marker for an empty volume; sparse ones fill over steps

### M-022 [minor, confirmed] `FindPointInCell` returns an error code as a cell index
- Tags: minor, confirmed, logic, quick
- Where: AVD.cpp:830-833 (commit 6124aa7d)
- What: `SETERRQ` inside a function returning `PetscInt` returns `PETSC_ERR_USER`, which is
  then used as I/J/K.
- Why it matters: a non-local marker silently corrupts data (`numMarkCell[]`, AVD.cpp:908)
  instead of stopping.
- Fix idea: return the index through an argument; return `PetscErrorCode`.
- Introduced by (git blame at HEAD):
  - Anton Popov — d3a5c7ee marker-surface-interaction major overhaul (2019-04-02) — added SETERRQ 'Non-local marker' inside PetscInt-returning function
  - Anton Popov — 5eae660a fucking bisection! (2015-11-19) — bisection search core of FindPointInCell
- At HEAD: still present, src/AVD.cpp:823-850 (FindPointInCell, SETERRQ at 830-833)
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-023 [minor, confirmed] RK2 deletes markers whose midpoint leaves the domain
- Tags: minor, confirmed, logic, quick
- Where: cvi.cpp:204-205 (commit 6124aa7d)
- What: `ADVelDeleteOutflow` runs on the midpoint position; a midpoint exactly at the upper
  wall also counts as outside (`x >= ex`).
- Why it matters: markers near walls with large dt are lost even if the final position is inside.
- Fix idea: clamp the midpoint to the domain.
- Introduced by (git blame at HEAD):
  - Adina Püsök — f848a86c CVI: stand-alone advection schemes and velocity interpolation routines (2015-07-20) — delete-outflow call on midpoint position; grank==-1 test (x>=ex outside)
- At HEAD: still present, src/cvi.cpp:204-205 (ADVelRungeKuttaStep) and ADVelDeleteOutflow 477-520
- Feedback requested from: Anton Popov
- Status: fixed (group C, low value)

### M-024 [minor, confirmed] `ADVGetSedPhase` uses stale cell numbers for passive tracers
- Tags: minor, confirmed, logic, quick
- Where: subgrid.cpp:673; LaMEMLib.cpp:731-734; passive_tracer.cpp:843 (commit 6124aa7d)
- What: called after `ADVExchange` but before `ADVMapMarkToCells`, so `actx->cellnum` no longer
  matches the markers.
- Why it matters: the sedimentation phase for passive tracers comes from wrong cells.
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 7968d4f4 Update with the old files and start working on it (2020-10-21) — calls ADVGetSedPhase in passive-tracer code using stale actx->cellnum
  - Andrea Piccolo — dd3b7cf1 I commit the first stage of the work (2020-10-21) — placed ADVAdvectPassiveTracer after ADVExchange with no ADVMapMarkToCells
- At HEAD: still present, src/LaMEMLib.cpp:718-721; src/passive_tracer.cpp:843; src/subgrid.cpp:631-ff. Note: ADVGetSedPhase itself (Pusok 521a14ec, 2015) is fine for the main markers; misuse is the passive-tracer call.
- Feedback requested from: Andrea Piccolo
- Status: open

### M-025 [minor, confirmed] Temperature-file setup loops over the wrong marker count
- Tags: minor, confirmed, memory, quick
- Where: marker.cpp:604-609, 620-642 (commit 6124aa7d)
- What: `ADVMarkSetTempFile` loops over `nmarkx*nmarky*nmarkz`, not `actx->nummark`. With
  `msetup = files` the counts can differ (write past nummark/markcap, or markers left unset).
  `Ix+1` reads out of range when xp == ex.
- Fix idea: loop to `actx->nummark`; clamp indices to [0, n-2].
- Introduced by (git blame at HEAD):
  - Arthur Bauville — 37bf8fbd Created ADVMarkSetTempFromFile function (2015-05-28) — wrote the loop and trilinear interpolation (Ix+1 unchecked)
  - Anton Popov — 4327aa93 ... (2017-01-18) — reworked nmarkx/y/z count lines and interpolation lines
- At HEAD: still present, src/marker.cpp:604-607, 614-642 (ADVMarkSetTempFile). Note: nummark = nmarkx*nmarky*nmarkz line itself is Popov 88b60a16 (2014).
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-026 [minor, confirmed] `ADVMarkSave` crashes when the save path has no '/'
- Tags: minor, confirmed, memory, quick
- Where: marker.cpp:294 (commit 6124aa7d)
- What: `(*strrchr(path, '/')) = '\0'` dereferences NULL for `mark_save_file = mdb`. The file
  is written with header -1, so APS is dropped. (Possible: restart with `msetup = files` may
  fail the exact-markers-per-cell check.)
- Introduced by (git blame at HEAD):
  - Anton Popov — 5bb51c5d test problem runs (2017-03-01) — added (*strrchr(path,'/'))='\0' without NULL check
- At HEAD: still present, src/marker.cpp:294
- Feedback requested from: Anton Popov
- Status: partly fixed (group A): APS still not saved

### M-027 [minor, confirmed] Periodic x-shift uses `>` while the rank test uses `>=`
- Tags: minor, confirmed, logic, quick
- Where: advect.cpp:1089 vs fdstag.h:327 (commit 6124aa7d)
- What: a marker with X[0] == ex is sent to the wrap neighbour but not shifted; that rank
  aborts with "Non-local point". Needs exact equality, so rare.
- Fix idea: `X[0] >= ex`.
- Introduced by (git blame at HEAD):
  - Anton Popov — 7b0ccd0e This belongs to previous commit (2022-02-06) — periodic wrap uses X[0] > ex
  - Anton Popov — 157753dd Bug fixes. (2017-03-24) — GET_POINT_REL_RANK uses x >= ex
- At HEAD: still present, src/advect.cpp:1089 vs src/fdstag.h:327. Note: Line 1089 re-blamed to 21a80d20 (2025-12-17 bugfix) only for a move; logic from 7b0ccd0e.
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-028 [minor, confirmed] Adiabatic gradient uses |z - Ztop|
- Tags: minor, confirmed, physics, quick
- Where: advect.cpp:2053; bc.cpp:2952 (commit 6124aa7d)
- What: `Z = PetscAbs(P->X[2]-Z_Top)`; rock markers above the reference level (topography)
  are heated as if they were below it.
- Fix idea: `Z = PetscMax(Z_Top - z, 0)`.
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 6cad036e Adiabatic gradient function (2021-01-26) — introduced PetscAbs(z - Z_Top) in both marker adiabat routines
- At HEAD: still present, src/advect.cpp:2053; src/bc.cpp:2952
- Feedback requested from: Andrea Piccolo
- Status: open

### M-029 [minor, confirmed] Ridge spreading rate can silently be 0
- Tags: minor, confirmed, input, quick
- Where: marker.cpp:1193-1200 (commit 6124aa7d)
- What: without `v_spread`, `bc->velin` is used, which is 0 without `bvel_velin` or with
  `velin_array`. The age becomes inf/huge, so the whole box gets T_top.
- Fix idea: require v_spread > 0.
- Introduced by (git blame at HEAD):
  - Jana Schierjott — 68d99ddd changed line 1099 to ridge->v_spread=PetscAbs(actx->jr->bc->velin) (2020-11-21) — fallback v_spread = |bc->velin| which can be 0
  - Jana Schierjott — e5fce5a2 added an additional input velocity for the ridge geometry (2021-03-16) — added optional v_spread (default 0), no check when both absent
- At HEAD: still present, src/marker.cpp:1193-1200
- Feedback requested from: Garrett Ito
- Status: open

### M-030 [minor, confirmed] Marker input not checked
- Tags: minor, confirmed, input, quick
- Where: parsing.cpp (`getIntParam` only checks max); marker.cpp:378, 536; advect.cpp:156-180;
  cvi.cpp:1239-1240 (commit 6124aa7d)
- What: negative phase IDs pass (out-of-range read at marker.cpp:536, caught later in
  `ADVCheckMarkPhases`); `nmark_sub <= 0`, `nmin > nmax`, `nmark_avd <= 0`, bad
  `nmark_lim` are accepted; `stagp_a = 0` is silently replaced by 2/3.
- Introduced by (git blame at HEAD):
  - Anton Popov — f17a7db6 Marker control parameters (2017-07-25) — nmark_lim/nmark_avd/nmark_sub read with no range checks
  - Tobias Baumann — f16da455 added temperature initialisation (linear gradient and phase based) (2017-11-21) — phase_temp[P->phase] indexing with unchecked (possibly negative) phase
  - Adina Püsök — f848a86c CVI: ... (2015-07-20) — stagp_a == 0 silently replaced by 2/3
- At HEAD: still present, src/advect.cpp:96-109,176-180; src/cvi.cpp:1238-1240; src/marker.cpp:536
- Feedback requested from: Anton Popov
- Status: fixed (group B, batch 1)

### M-031 [minor, confirmed] Memory leak in `ADVelReAllocStorage`
- Tags: minor, confirmed, memory, quick
- Where: cvi.cpp:800, 826-828 (commit 6124aa7d)
- What: `vi->cellnum` is freed but `vi->markind` is overwritten without `PetscFree`.
- Introduced by (git blame at HEAD):
  - Adina Püsök — f848a86c CVI: stand-alone advection schemes and velocity interpolation routines (2015-07-20) — allocated vi->markind in realloc without freeing the old array
- At HEAD: still present, src/cvi.cpp:800,825-828 (ADVelReAllocStorage). Note: Blame via -C points to 521a14ec (Pusok, 2015-05-28) as the copied pattern origin.
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-032 [minor, confirmed] Dikes without free surface read `phRat[-1]`
- Tags: minor, confirmed, memory, quick
- Where: src/dike.cpp:352 (`GetDikeContr`), 673 (`Compute_sxx_magP`); surf.cpp:37
  (commit 6124aa7d)
- What: `AirPhase` is -1 when `surf_use = 0`, and is used as an index without a check.
- Why it matters: out-of-bounds read; `dikeRHS` and the stress integral get a random extra
  term.
- Fix idea: guard with `AirPhase >= 0`, or require a free surface when `act_dike = 1`.
- Introduced by (git blame at HEAD):
  - Garrett Ito — 20b0a12f dike.cpp: if cell contains air and dike, then give full divergence (2022-03-26) — phRat[i]+phRat[AirPhase] with AirPhase possibly -1
  - Garrett Ito — dbf716ed dike.cpp: Pattern of sxx_ave and dPm vs x match Paraview output (2022-07-31) — svCell->phRat[AirPhase] < 1.0 with unchecked AirPhase
- At HEAD: still present, src/dike.cpp:352 (GetDikeContr), 673 (Compute_sxx_magP); src/surf.cpp:37. Note: Line 352 blames to Jana 335839e0 (comment-only); logic from Ito 20b0a12f. AirPhase=-1 default is Popov 94177af5 (2017).
- Feedback requested from: Garrett Ito
- Status: open

### M-033 [minor, possible] Dynamic dike can produce NaN
- Tags: minor, possible, logic, quick
- Where: src/dike.cpp:1095, 753-754 (commit 6124aa7d)
- What: the local dike azimuth is `atan(sumslope/sumadd)`. If no other dike cell lies
  within `filty` (cell size in y larger than `filty`, default 1.5 km; or a one-row dike),
  `sumadd = 0` and `azim = NaN`, which spreads to all smoothed stresses. In
  `Compute_sxx_magP`, a column with no rock below `Tsol` has `liththick = 0` → divide by 0.
- Fix idea: fall back to `azim = 0` when `sumadd == 0`; guard `liththick > 0`.
- Introduced by (git blame at HEAD):
  - Garrett Ito — c2f25559 Updates since last sync with bitbucket/lamem kaus (2023-08-06) — azim = atan(sumslope/sumadd) without sumadd==0 guard; divides by liththick
  - Garrett Ito — de4c4637 dike.cpp: now using Petsc..Scalar math tools (2025-11-07) — only renamed atan to PetscAtanScalar
- At HEAD: still present, src/dike.cpp:1095 (azim), 753-754 (liththick division). Note: de4c4637 is cosmetic; real introducer is c2f25559.
- Feedback requested from: Garrett Ito
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
- Introduced by (git blame at HEAD):
  - Garrett Ito — cf83ee4d dike.cpp: created DynamicDike_ReadRestart, DynamicDike_WriteRestart (2023-03-18) — added TSSolCreate(ts, fb) in restart reader, resetting dt
  - Garrett Ito — 2847e43e Smooth_Sxx_eff does time averaging over istep_nave steps (2023-01-13) — DBDikeCreate sets istep_count=istep_nave, reset on restart
- At HEAD: still present, src/dike.cpp:1588 (TSSolCreate in DynamicDike_ReadRestart), 219 (istep_count reset)
- Feedback requested from: Garrett Ito
- Status: open

### M-035 [minor, confirmed] Dike opening rate and units
- Tags: minor, confirmed, input, quick
- Where: src/dike.cpp:306-345 (and 414-453), 212, 855, 1230 (commit 6124aa7d)
- What: dike divergence uses `|bc->velin|` as spreading rate. Without `bvel_velin`
  (e.g. spreading set by other BCs) it is 0 and the dike silently does nothing.
  `magPwidth` is not divided by `scal->length` (`filtx`, `filty` are).
- Fix idea: make the spreading rate an explicit dike input or check it is > 0; scale
  `magPwidth`.
- Introduced by (git blame at HEAD):
  - Jana Schierjott — 6317b311 moved the function GetDikeContr to dike.cpp (2021-06-30) — spreading rate taken from |bc->velin| (idea since 3b1fb801, 2020-12-03)
  - Garrett Ito — c2f25559 Updates since last sync with bitbucket/lamem kaus (2023-08-06) — added magPwidth read with scale 1.0 (not length)
- At HEAD: still present, src/dike.cpp:306,317,340,414,425,448 (v_spread=|bc->velin|); 212,1230 (magPwidth). Note: Velin-based spreading originated in Jana's dike prototype in constEq/JacRes.
- Feedback requested from: Garrett Ito
- Status: open

### M-036 [minor, confirmed] Tracer melt phase search checks the wrong marker
- Tags: minor, confirmed, logic, memory, quick
- Where: src/passive_tracer.cpp:599, also 370, 607, 1049 (commit 6124aa7d)
- What: in the loop over markers of the cell, the phase test uses `actx->markers[ii].phase`
  (loop counter) instead of `actx->markers[id_m].phase`. If no marker passes the test, `dist` is
  empty and `dist.begin()->second` reads past the end. The same empty-list read happens in
  `ADV_Assign_Phase` and `Check_advection_condition` if a cell has no markers.
- Fix idea: use `id_m`; check `dist.size()` before use.
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 14e0f271 Update with the old files and start working on it (2020-10-21) — mat[actx->markers[ii].phase] uses loop counter; unguarded dist.begin()
  - Yun Luo — e4778c61 New coupling module based on LaMEM 3.0.1 (2026-08-02) — touched line 599 (pdn test style) but kept markers[ii]
- At HEAD: still present, src/passive_tracer.cpp:599 (ii vs id_m), 370, 606-607, 936, 1048-1049. Note: Original bug line in 2020 Piccolo commits (also beeec8c5 2020-11-10); Yun's change was cosmetic.
- Feedback requested from: Andrea Piccolo, Yun Luo
- Status: open

### M-037 [minor, confirmed] CFL and CFLMAX range checks can never trigger
- Tags: minor, confirmed, input, quick
- Where: src/tssolve.cpp:56, 61 (commit 6124aa7d)
- What: `if(ts->CFL < 0.0 && ts->CFL > 1.0)` is never true; should be `||`. Same for CFLMAX.
  CFL > 1 or negative values are accepted.
- Introduced by (git blame at HEAD):
  - Anton Popov — 5b941b56 Hexahedron geometric primitives for phase assignment (2017-03-07) — wrote CFL/CFLMAX range checks with && instead of ||
- At HEAD: still present, src/tssolve.cpp:56,61
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-038 [minor, confirmed] Automatic coarse-grid CPU count uses the local cell count
- Tags: minor, confirmed, parallel, quick
- Where: src/options.cpp:482 (commit 6124aa7d); fdstag.cpp:1645 (local sizes)
- What: with `coarse_num_cpu = -1` the target is `ceil(coarse_num_local_cells / coarse_cells_per_cpu)`.
  `coarse_num_local_cells` is the count on one rank, so the result is almost always 1 and the
  whole coarse grid goes to one CPU, however large it is.
- Fix idea: multiply by the number of ranks (total coarse cells).
- Introduced by (git blame at HEAD):
  - Anton Popov — c9dcd497 ... (2025-08-25) — changed coarse CPU estimate from global cell count to coarse_num_local_cells
- At HEAD: still present, src/options.cpp:482 (+ fdstag.cpp:1666-1671 local size). Note: Earlier b7f5b864 version used global Nx*Ny*Nz; local count regression introduced in c9dcd497.
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Boris Kaus — eeeee733 the passive_tracer output files were not being written. It requires an (2020-11-20) — added out_ptr read directly into ctrl.Passive_Tracer
  - Andrea Piccolo — 8fe9dd92 paraview passive tracers (2020-11-18) — original Passive_Tracer-based PVPtr output switch that out_ptr overwrote
- At HEAD: still present, src/paraViewOutPassiveTracers.cpp:32. Note: Passive_Tracer flag itself and PVPtrCreate early-return come from Andrea Piccolo (8fe9dd92, 2020-11-18).
- Feedback requested from: Andrea Piccolo, Boris Kaus
- Status: open

### M-040 [minor, confirmed] `UpdatePVDFile` uses the file name after `free` in its error message
- Tags: minor, confirmed, memory, quick
- Where: src/paraViewOutBin.cpp:770-772 (commit 6124aa7d)
- What: `free(fname)` comes before `if(fp == NULL) SETERRQ(..., "cannot open file %s", fname)`.
- Why it matters: when the .pvd file can't be opened (e.g. restart in a folder without the old
  .pvd file, which opens it with "r+b"), the error message reads freed memory, which can print
  garbage or crash.
- Fix idea: move `free(fname)` after the check (as in the other open calls).
- Introduced by (git blame at HEAD):
  - Anton Popov — f1b6ede1 This belongs to previous commit (2017-03-15) — moved free(fname) before the fp==NULL error using fname
- At HEAD: still present, src/paraViewOutBin.cpp:770-772
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Anton Popov — f1b6ede1 This belongs to previous commit (2017-03-15); 2f4b267b Binary specificator for file opening (2022-09-09) — restart-mode r+b open and offset seek without truncating the .pvd
- At HEAD: still present, src/paraViewOutBin.cpp:766-797. Note: Missing-truncate bug; blamed code is the r+b/ttime-based open and offset logic.
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### M-042 [minor, confirmed] `Adjoint_ScaleCostFunction = Var` uses only the last point
- Tags: minor, confirmed, logic, quick
- Where: src/adjoint.cpp:920-947
- What: `var = pow(...)` instead of `+=`. The mean uses `abs(sum)`, so the variance is taken
  around the wrong value for negative data. `mdI = 0` (Ap = 2/3) or zero mean divides by zero.
  "Mean" divides the squared misfit by the mean (not its square) in non-dimensional units, so
  the result still carries units.
- Introduced by (git blame at HEAD):
  - Georg Reuber — e1d9b069 Cost function evaluation is a separate function; introduced cost funct (2020-05-22) — wrote var = pow(...) (no +=) and abs(mean) in Var scaling
- At HEAD: still present, src/adjoint.cpp:937-947
- Feedback requested from: Boris Kaus
- Status: open
- Boris: seems a bug; fix

### M-043 [minor, confirmed] Scaling-law exponents are wrong for `log10` parameters
- Tags: minor, confirmed, physics, quick
- Where: src/adjoint.cpp:3024-3040, 3080-3095, 3151
- What: `b = grad * P / F` assumes `grad = dF/dP`. For `log10 = 1`, `P` is log10(p) and `grad`
  is `dF/dlog10(p)`, so `b` should be `grad / (F ln 10)`, and `A` should use `p`, not `log10(p)`.
  "Velocity check" is always equal to F by construction. The file writes non-dimensional
  coordinates.
- Introduced by (git blame at HEAD):
  - Boris Kaus — ec790532 Added printing of scaling laws, which can be activated with: (2020-05-24); da491451 added the possibility to do inversion with the log10(parameter), which (2020-05-31) — b = grad*P/F scaling law; log10 parameters added later without adapting it
- At HEAD: still present, src/adjoint.cpp:3026,3043,3109 (ScalingLaws)
- Feedback requested from: Boris Kaus
- Status: open
- Boris please fix

### M-044 [minor, confirmed] A parameter value of 0 gives a zero FD step
- Tags: minor, confirmed, logic, quick
- Where: src/adjoint.cpp:1866 (brute-force FD), 1995 (adjoint dr/dp), 1328 (descent start step)
- What: steps are relative (`eps * value`). A parameter equal to 0 (e.g. a `log10` value of 0,
  or cohesion 0) divides by zero. A zero gradient gives an infinite start step in the built-in
  descent.
- Fix idea: use `eps * max(|value|, 1)` or an absolute step.
- Introduced by (git blame at HEAD):
  - Boris Kaus — 5d65e2ec - Scaling laws and exponents cannot be computed if the value of a para (2020-05-24); 8e93a88f WIP: attempt to compute dr/dp by changing command-line options. (2020-05-15) — relative FD steps Perturb = CurVal*eps (zero for zero value)
  - Georg Reuber — 3492efcd fixing a bunch of things in the GD algorithm (2020-06-08) — start step factor2array = Scale_Grad*Par/|grad|
- At HEAD: still present, src/adjoint.cpp:1328,1866,1995
- Feedback requested from: Boris Kaus
- Status: open
- Boris: fix

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
- Introduced by (git blame at HEAD):
  - Jana Schierjott — 85e29cfc removed the dikeRHS field from being created for paraview because it i (2021-06-03); 12491b0f changes for strain rate removal in JacREs.cpp and removed the old vers (2021-04-28) — in-place dike strain-rate removal on owned cells only
- At HEAD: still present, src/JacRes.cpp:1008-1016
- Feedback requested from: Garrett Ito
- Status: open

### M-046 [minor, possible] Dike opening still enters stabilization stress and shear heating
- Tags: minor, possible, physics, medium
- Where: src/constEq.cpp:798-822 (called from JacRes.cpp:1079)
- What: cellConstEq uses `svCell->dxx`, which is not dike-corrected, for the stabilization
  stress `2*eta_st*dxx` and for the dissipative strain rate in `Hr`. The corrected value is
  used only for the main stress.
- Why it matters: dike opening creates false stabilization stress and false shear heating
  in dike cells.
- Introduced by (git blame at HEAD):
  - Jana Schierjott — 12491b0f changes for strain rate removal in JacREs.cpp and removed the old vers (2021-04-28); 85e29cfc removed the dikeRHS field from being created for paraview because it i (2021-06-03) — dike-corrected dxx passed separately; svCell->dxx left uncorrected
  - Anton Popov — 49845492 Stabilizations and regularizations. (2020-05-10); 8890a329 Update (2020-03-03) — stabilization stress and Hr use svCell->dxx
- At HEAD: still present, src/constEq.cpp:799-801,817-824; JacRes.cpp:1011-1016. Note: Bug arises from dike correction (Jana) interacting with older Popov stabilization/heating code.
- Feedback requested from: Anton Popov, Garrett Ito
- Status: open

### M-047 [minor, possible] `init_lith_pres` only sets the pressure history; weak stop test
- Tags: minor, possible, logic, medium
- Where: src/JacRes.cpp:1545-1703
- What: the computed lithostatic pressure is copied to `svBulk.pn` and to the markers, but
  not to the pressure in `gsol`, which stays 0. The convergence test compares global sums
  of |p|. Local changes can cancel out. If the sum is 0 (no gravity), you get 0/0 = NaN,
  it runs all iterations and prints a false warning.
- Why it matters: the option name suggests the solution pressure is initialized, but it
  isn't. Check whether that's intended.
- Introduced by (git blame at HEAD):
  - Anton Popov — 84717c72 fixes in iterative lithostatic pressure (2019-08-22); 8afca1f8 Update (2020-03-02) — norm-sum stop test; copies p only to svBulk.pn and markers
  - Arne Spang — c87cec03 Added initialization of pressure with lithostatic pressure to stabiliz (2019-08-22) — original init_lith_pres feature
- At HEAD: still present, src/JacRes.cpp:1545-1700
- Feedback requested from: Anton Popov
- Status: open

### M-048 [minor, confirmed] Phase temperature `T = 0` (°C) is treated as "not set"
- Tags: minor, confirmed, input, medium
- Where: src/phase.cpp:798-803; src/marker.cpp:523-536
- What: `T`, `T_Nu`, `T_liq` and `T_sol` are shifted to Kelvin only when they are non-zero,
  and 0 means "unset". A user who sets air temperature `T = 0` [°C], which is common, gets
  no phase temperature at all. The code already warns for the same problem with `TRef_fk`
  (phase.cpp:703).
- Fix idea: use a DBL_MAX sentinel for "unset", like `gwLevel`.
- Introduced by (git blame at HEAD):
  - Anton Popov — 2d91a3bf Temperature initialization (started) (2017-11-20); cff764ed ps-dc creep update (should be functional for quasi-harmonic mean) (2019-06-25) — T!=0 treated as set; Kelvin shift only if nonzero
  - Tobias Baumann — f16da455 - added temperature initialisation (linear gradient and phase based) - (2017-11-21) — phase-based marker T assignment skipping zero
  - Garrett Ito — f1da6107 Wrote code for Behn & Ito Dike heat terms but haven't compiled phase.c (2021-07-27) — same non-zero shift for T_liq/T_sol
- At HEAD: still present, src/phase.cpp:786-791; src/marker.cpp:523-536
- Feedback requested from: Anton Popov, Garrett Ito
- Status: open

### M-049 [minor, possible] Moduli averaged with non-elastic phases counted as 0
- Tags: minor, possible, physics, medium
- Where: src/constEq.cpp:603-626 (G), :697 and :756 (Kb)
- What: `Gavg = sum(phRat*G)` counts phases without G as 0. In a cell that's 50% elastic rock
  and 50% non-elastic phase (often air), I2Gdt doubles. The phase viscosity, though, uses
  each phase's own G. Kb works the same way: a mixed cell becomes more compressible than
  either phase.
- Why it matters: elastic stress history is too strong in mixed cells, for example at the
  free surface. It's hard to see, but it isn't consistent.
- Introduced by (git blame at HEAD):
  - Anton Popov — a4b35853 Ported recent stuff from master (2017-02-13) — Gavg/Kavg summed with non-elastic phases counted as zero
- At HEAD: still present, src/constEq.cpp:620-623 (G), 697,756 (Kb)
- Feedback requested from: Anton Popov
- Status: open

### M-050 [minor, confirmed] Slope erosion and topographic diffusion ignore periodic x
- Tags: minor, confirmed, logic, medium
- Where: src/surf.cpp:1167, 1223, 1643
- What: `I1 < 0 -> I1 = I` also runs for periodic grids, so the left edge is one-sided (no
  x-diffusion there), while the right edge wraps through the ghost node. FreeSurfAdvectTopo
  handles `periodic` correctly; these two functions don't.
- Introduced by (git blame at HEAD):
  - Wenrong Cao — 1a44f790 Add topographic diffusion feature to free surface module (2026-04-04) — I1<0 -> I1=I clamp without periodic check (slope erosion, topo diffusion)
- At HEAD: still present, src/surf.cpp:1167,1223,1643
- Feedback requested from: Wenrong Cao
- Status: open

### M-051 [minor, confirmed] Rank-0-only errors can hang; restart reads aren't checked
- Tags: minor, confirmed, parallel, restart, medium
- Where: src/tools.cpp:293-296 (DirMake), 319-322, 342-345; tools.cpp:80, fdstag.cpp:292-294
- What: DirMake/DirRemove/DirRename raise SETERRQ on rank 0 only. In DirMake the other ranks
  are waiting in MPI_Barrier, so a failed mkdir hangs the run instead of stopping it. Restart
  readers ignore the return value of fread, so a truncated restart file loads garbage without
  any message.
- Fix idea: broadcast the status and fail on all ranks; check fread counts.
- Introduced by (git blame at HEAD):
  - Anton Popov — 1a6c9c56 # (2017-02-20); cd2c66d8 Upgrade to petsc-3.18 (2023-02-28) — rank-0-only SETERRQ in DirMake/DirRemove/DirRename
  - Anton Popov — e918c7d4 JacRes (2016-11-23); 1d3fae78 no comments (2015-07-27) — unchecked fread in restart readers
- At HEAD: still present, src/tools.cpp:295,321,344,80; src/fdstag.cpp:292-294
- Feedback requested from: Anton Popov
- Status: fixed (batch 4): no hang, PETSc aborts all ranks; fread counts checked

### M-052 [minor, confirmed] Sign of `velin` on Right/Back faces flips depending on `bvel_face_out`
- Tags: minor, confirmed, input, medium
- Where: src/bc.cpp:2103-2123 (x), 2167-2187 (y)
- What: with `bvel_face_out` set, Right/Back use `-vel` (positive `velin` = inflow). Without
  it, they use `+vel` (positive `velin` = outflow). The mass-balance `velout` uses the same
  formula in both cases.
- Why it matters: the same input means opposite flow directions; easy to set up wrong.
- Fix idea: use one convention (inflow positive) in both branches, and document it.
- Introduced by (git blame at HEAD):
  - Iskander Ibragimov — d6e7958a Added new inflow from boundary condition options (2023-05-17) — added face_out branch with -vel on Right/Back faces
  - Andrea Piccolo — 1bc3ab4d Testing the new boundary conditions, still working on it (2020-10-04) — non-face_out branch (+vel on Right/Back faces) from earlier BC code
- At HEAD: still present, src/bc.cpp:2105-2114 (x), 2169-2179 (y) vs 2122-2123, 2186-2187
- Feedback requested from: Andrea Piccolo
- Status: open

### M-053 [minor, confirmed] Fixed-cell file read: unchecked I/O, one-rank errors hang
- Tags: minor, confirmed, parallel, input, medium
- Where: src/bc.cpp:1229-1274 (BCReadFixCell), 1148 (BCReadRestart)
- What: each rank opens its own file; `SETERRQ(PETSC_COMM_WORLD, ...)` on one rank (missing
  file, wrong size) can hang the others. `stat` and `fread` results are not checked; `fp`
  is not closed on the size error. Same pattern as M-051.
- Fix idea: collect the status with an MPI reduction and fail on all ranks; check counts.
- Introduced by (git blame at HEAD):
  - Anton Popov — 139c4a1a constrained cells form parallel files (2017-08-02) — per-rank fixed cell file read; stat/fread unchecked
- At HEAD: still present, src/bc.cpp:1148, 1257-1269
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Thomas Morrow — a6a42c90 Added matching outflow boundary conditions on side with compensating i (2020-10-30); 5a7ecfab Listed variables for balanced inflow are now OPTIONAL. A test model co (2020-11-01) — compensating inflow face 5 without mass-balance check
  - Panupong Kongpet — 7d0eb7a1 Add velin_net_num_periods, velin_net_time_delims, and bvel_velin_net f (2026-01-25) — added velin_net shift on right side
- At HEAD: still present, src/bc.cpp:2129-2142, 2213-2223
- Feedback requested from: Anton Popov (assigned by code area; none of the main suspects wrote these lines)
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
- Introduced by (git blame at HEAD):
  - Anton Popov — 94ea1409 Tangential boundary velocity support (2024-05-07); eaa57323 Assemble diagonal directly into vector (matrix-free) (2025-05-05) — matrix-free two-point constraints and diagonal flags at i==0/mnx ignore periodic
  - Anton Popov — 47bbc92f wBFBT support for periodic BC (completed) (2025-10-31) — wBFBT periodic handling leaves cfv[0],cfv[1] at 1
- At HEAD: still present, src/matFree.cpp:639-640,683-684,1307-1308,1358-1359; src/matBFBT.cpp:208-209
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Anton Popov — 0b793ef0 update ... (2025-04-20) — coarse BC vectors set to DBL_MAX, ghost constraints not restricted
- At HEAD: still present, src/matData.cpp:642-645, 808-811. Note: Blame shows later PetscCall cleanup 1856f403; originating commit via git log -S.
- Feedback requested from: Anton Popov
- Status: partly fixed (group B, batch 5): 3D fixed, 2D left (needs matrix-free MG)

### M-057 [minor, confirmed] ADVCheckCorners depends on CPU count; same seed; centre-plane gap
- Tags: minor, confirmed, parallel, medium
- Where: advect.cpp:1368-1376, 1412-1413, 1459, 1568 (commit 6124aa7d)
- What: the parent is searched in local cells only; a new `PetscRandom` with the default seed
  is made every call on every rank; a marker with xp == xc counts in no corner, so markers
  are injected into non-empty corners; the count is printed for rank 0 only.
- Fix idea: use ghost markers or document; seed once per rank; use `>=` on one side; reduce
  the count.
- Introduced by (git blame at HEAD):
  - Adina Püsök — 766f54c3 1. Marker control (solves the "Empty Control Volume" error) * activate (2015-03-30); e5016a08 added local neighbor cells influence for corner marker injection to de (2015-06-15) — ADVCheckCorners: local-cell search, strict corner tests, per-call random context
  - Anton Popov — d3a5c7ee marker-surface-interaction major overhaul (2019-04-02) — switched random context to PETSC_COMM_SELF
- At HEAD: still present, src/advect.cpp:1368-1376,1412-1413,1459,1568
- Feedback requested from: Anton Popov
- Status: fixed (group C, low value)

### M-058 [minor, confirmed] AVD half-volume rule uses processor boundaries
- Tags: minor, confirmed, parallel, medium
- Where: AVD.cpp:747-749, 791-793 (commit 6124aa7d)
- What: `(i == 0) | (i+1 == mv->M)` uses the local count, so a face volume cut by a processor
  boundary is checked as two halves needing nmin/2+1 each. `nmax` is not halved for half
  volumes.
- Why it matters: marker control depends on CPU count; boundary volumes can hold twice the
  density before deletion.
- Introduced by (git blame at HEAD):
  - Anton Popov — 8763ccaf Multiple updates from concurrent branches (2017-01-11) — half-volume rule uses local mv->M/N/P index bounds
- At HEAD: still present, src/AVD.cpp:747-749,791-793. Note: Introducing commit is a large multi-branch update by Popov (new AVD marker control).
- Feedback requested from: Anton Popov
- Status: fixed (AVD batch)

### M-059 [minor, confirmed] AVD "edge" control volumes are face volumes; labels swapped
- Tags: minor, confirmed, logic, medium
- Where: AVD.cpp:944-946, 865-884, 813-816 (commit 6124aa7d)
- What: for each dir only one coordinate is shifted to cell centres, so the volume is centred
  on a face, not an edge. The log prints dir 0 as "XYED" and dir 2 as "YZED", the reverse of
  lines 695-697.
- Why it matters: volumes around shear-stress points are never checked; log is misleading.
- Introduced by (git blame at HEAD):
  - Anton Popov — 8763ccaf Multiple updates from concurrent branches (2017-01-11) — face-centred volumes called edges; XYED/YZED labels swapped
- At HEAD: still present, src/AVD.cpp:944-946,865-870,814-816 vs 695-697
- Feedback requested from: Anton Popov
- Status: fixed (AVD batch)

### M-060 [minor, confirmed] Polygon file sizes not checked against buffers
- Tags: minor, confirmed, memory, medium
- Where: marker.cpp:1466, 1622-1624, 1664-1744 (commit 6124aa7d)
- What: `Polys[_max_polygons_]` (10) overflows with more polygons in one slice; `Vol.num`,
  `len`, `Fcount` are not checked against the allocated sizes; VLAs `SxAll[Vol.num]`.
- Why it matters: a bad or old geomIO file overwrites memory.
- Introduced by (git blame at HEAD):
  - Arne Spang — 82be9435 Extended ability to read geomIO polygons (2022-05-05) — polygon file read without size checks against Polys[]/VLAs
- At HEAD: still present, src/marker.cpp:1466,1643,1676-1677,1730-1732
- Feedback requested from: Anton Popov
- Status: fixed (group B, batch 2)

### M-061 [minor, confirmed] Corner output at the periodic x seam uses one side only
- Tags: minor, confirmed, logic, medium
- Where: interpolate.cpp:102, 158, 220, 377 (commit 6124aa7d)
- What: in periodic mode `i == mx` never happens, so only the `i == 0` clamp fires and the
  corner value at x = 0 ignores ghost cell -1. Output only; velocity (`use_bound = 1`) is fine.
- Introduced by (git blame at HEAD):
  - Anton Popov — 88b60a16 code integration & functionality restore update version 0.0 (2014-09-09) — corner interpolation clamps only at i==0 / i==mx, no periodic handling
- At HEAD: still present, src/interpolate.cpp:102-103,158-159,220-222,377-378. Note: Original 2014 code; periodic x added later without updating this.
- Feedback requested from: Anton Popov
- Status: fixed (x clamps skipped when periodic; ghost -1 holds the wrapped value)

### M-062 [minor, confirmed] Halfspace temperatures use hard-coded kappa = 1e-6
- Tags: minor, confirmed, physics, medium
- Where: marker.cpp:1013, 1152, 1215; phase_transition.cpp:1197, 1339 (Box halfspace T)
  (commit 6124aa7d)
- What: ignores phase conductivity and heat capacity; in non-dimensional mode it becomes a
  non-dimensional 1e-6.
- Fix idea: make kappa an input, or reject in non-dimensional setups.
- Introduced by (git blame at HEAD):
  - Boris Kaus — d92a75d9 added the possibility to add an initial temperature field or thermal s (2019-01-17); d7142fe3 Finished implemention of Box-like phase transition, which also allows  (2020-11-17) — hard-coded kappa=1e-6 for halfspace T (layer, box, PT box)
  - Jana Schierjott — 993280d8 option added for oblique and perpendicular ridges. SetPhaseRidge, Comp (2020-11-17) — same hard-coded kappa for ridge
- At HEAD: still present, src/marker.cpp:1014,1152,1215; src/phase_transition.cpp:1197,1339
- Feedback requested from: Boris Kaus, Garrett Ito
- Status: open
- Boris:yes,  please fix (easy)

### M-063 [minor, confirmed] Rank-local hard errors can hang; sedimentation-phase error on fast uplift
- Tags: minor, confirmed, parallel, medium
- Where: subgrid.cpp:85, 606-608; advect.cpp:1906; tools.cpp:293, 798 (commit 6124aa7d)
- What: "Empty control volume", "wrong phase" and `DirMake` errors are raised on one rank in
  collective code, so other ranks hang. subgrid.cpp:606 aborts with "Incorrect sedimentation
  phase" when the surface rises more than one cell per step or has overhangs.
- Introduced by (git blame at HEAD):
  - Anton Popov — 250306d9 subgrid update (2019-03-14); 60dc961e marker-free-surface-interaction update (2019-03-26); d3a5c7ee marker-surface-interaction major overhaul (2019-04-02) — subgrid errors raised on one rank; sedimentation phase error
  - Adina Püsök — 848cbfc3 added error check for marker phases (2015-06-17) — wrong-phase SETERRQ_SELF on single rank
  - Anton Popov — 1a6c9c56 # (2017-02-20) — DirMake rank-0 only error
- At HEAD: still present, src/subgrid.cpp:85,606-608; src/advect.cpp:1906-1908; src/tools.cpp:295,798-800
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Boris Kaus — 4a64c30d added the possibility to set plastic strain to zero as part of the pha (2020-12-01); eb789229 bugfix (2021-12-16) — APS reset every step on reset side; PhaseOutside=-1 branch
  - Jana Schierjott — af231bf1 added function in Transition() to be able to not turn air particles ba (2021-04-01); d7ac9cd6 additional changes for the new phase transition box that excludes the  (2021-06-08) — NotInAirBox check never sets InAbove
- At HEAD: still present, src/phase_transition.cpp:886-905,922,1024-1029
- Feedback requested from: Boris Kaus, Garrett Ito
- Status: open
- Boris: yes, fix

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
- Introduced by (git blame at HEAD):
  - Boris Kaus — f3527d46 PhaseTransitions: added a possibility to only check particles nearby a (2021-12-13) — BoxVicinity with dzBox computed from y bounds
  - Garrett Ito — 1a09a4df phase_transition:  compiled, but not producing dike phase on Mahi, and (2022-06-20) — bounds[] only set when nsegs not given
- At HEAD: still present, src/phase_transition.cpp:851-864 (BoxVicinity), 398-412 (bounds)
- Feedback requested from: Garrett Ito, Boris Kaus
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
- Introduced by (git blame at HEAD):
  - Arne Spang — 1c27c92e New options that allow planning the timesteps for the whole simulation (2022-05-04) — schedule built without limit check; nstep_max cut; shift loop off by one
- At HEAD: still present, src/tssolve.cpp:407-419,457-462
- Feedback requested from: Anton Popov
- Status: fixed (group B, batch 2)

### M-067 [minor, confirmed] Time-step restart repeats phase transitions
- Tags: minor, confirmed, logic, medium
- Where: src/LaMEMLib.cpp:713 (`if(restart) continue;`), src/phase_transition.cpp:967-968 (commit 6124aa7d)
- What: when a fixed time step exceeds CFLMAX, the loop starts the step again from the top.
  `Phase_Transition` runs again, so `MovingBox` moves the box a second time (old dt, then new dt),
  and dike time averaging counts the step twice. With adjoint modes, the gradient for the
  rejected step was already computed and added.
- Fix idea: do the CFL check before applying step-start changes, or skip them on a restart.
- Introduced by (git blame at HEAD):
  - Jana Schierjott — ddecc04c added the new phase transition Box moving routine, including tests (2021-10-26) — MovingBox moves box each time step call
  - Garrett Ito — deae58d8 MovingBox:  Now moves PhaseTransition->celly_xboundL, xboundR.  Passes (2022-06-26) — MovingBox moves celly_xboundL/R by v_box*dt
  - Anton Popov — 03950698 @ (2017-02-22); 9a108627 Explicit time stepping + elasticity + kinematic bc block. Partially po (2017-03-14) — CFL restart "continue" after step-start transitions
- At HEAD: still present, src/LaMEMLib.cpp:697-700; src/phase_transition.cpp:958-970. Note: Interaction bug: Popov CFL restart loop (2017) vs later Jana/Garrett stateful MovingBox; first author arguable.
- Feedback requested from: Anton Popov, Garrett Ito
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
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 14e0f271 Update with the old files and start working on it (2020-10-21); 00dc1a95 Big changes: i try to simplify the function structure, introducing all (2020-11-11) — ID formula with ny/nx swapped; silent default Always; ActiveValue==0 rejected
  - Mitchell McMillan — 924fda26 add optional APS output to passive tracers; defaults to 0 (2025-08-27) — APS copied once from closest marker, never advected
- At HEAD: still present, src/passive_tracer.cpp:262,376,60-63,523 (aps unused)
- Feedback requested from: Andrea Piccolo
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
- Introduced by (git blame at HEAD):
  - Anton Popov — 94177af5 . (2017-02-16); d8d709fa Direction vectors for:  * Maximum horizontal compression stress (SHmax (2015-12-01) — x>=0 sign rule on axial directions; averaged linearly to corners
- At HEAD: still present, src/JacResAux.cpp:77,156; src/outFunct.cpp:742-743,767-768
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Boris Kaus — e285ed37 added the possibility to select all material parameters within the Par (2020-05-23); 1748467f Added the possibility to exclude particular phases if the option "Type (2020-06-15); 0fcb16a9 - added more orinting statements, to clarify which parameters are comp (2020-05-29) — scan loop without blockID++, unchecked ExcludePhase/array limit, par_str leak
- At HEAD: still present, src/adjoint.cpp:204-260,296-400
- Feedback requested from: Boris Kaus
- Status: open
- Boris; fix

### M-071 [minor, confirmed] Adjoint solve leaves the forward KSP on the `as_` prefix and cells in a perturbed state
- Tags: minor, confirmed, logic, medium
- Where: src/adjoint.cpp:1935, 1947, 1995-2024
- What: the adjoint uses the SNES's own KSP and sets its prefix to `as_`. It is never set back,
  and the adjoint runs every time step, so all later forward solves read `as_` options. After the
  parameter loop the material DB is restored, but `svCell/svEdge` still hold values from the
  last perturbed residual (`FormResidual`); these go to output and the next step (tiny effect,
  eps = 1e-6).
- Fix idea: use a separate KSP for the adjoint; recompute the residual after the loop.
- Introduced by (git blame at HEAD):
  - Boris Kaus — cc0c2404 fixed the remaining memory leak in gradientdescent (2020-05-18); 8e93a88f WIP: attempt to compute dr/dp by changing command-line options. (2020-05-15) — adjoint reuses SNES KSP with as_ prefix; perturbed FormResidual leaves cell state
  - Georg Reuber — 67738f08 Added PSD angle derivative calculation (2020-06-02) — second as_ prefix branch (MfitType 1)
- At HEAD: still present, src/adjoint.cpp:1934-1935,1946-1947,1995-2024
- Feedback requested from: Boris Kaus
- Status: open
- Bois: try fixing, low priority

### M-072 [minor, confirmed] `SyntheticForwardRun` saves a zero vector; the `OFdef = 0` path doesn't exist
- Tags: minor, confirmed, cleanup, medium
- Where: src/adjoint.cpp:1123-1132, 62-66, 2742; LaMEMLib.cpp:690
- What: in `SyntheticForwardRun` the adjoint objective is never called, so `IOparam->xini` stays
  zero and `Forward_Solution_Vel.bin` is all zeros (it would only hold point targets anyway).
  The header describes loading `Forward_Solution.bin` with `OFdef = 0`, but nothing loads it,
  and `Adjoint_ObjectiveFunctionDef` can only be 1.
- Fix idea: save `jr->gsol` and add the load, or remove the mode and the comment.
- Introduced by (git blame at HEAD):
  - Georg Reuber — cf0ea953 Added the adjoint gradient computation and the TAO framework to LaMEM  (2018-01-11); fdbed43e Clean up of LaMEM.cpp. (2018-01-15); 2da29341 Initial creation of the branch with first working version (2020-05-08) — SyntheticForwardRun saves xini; OFdef=0 documented but not implemented
- At HEAD: still present, src/adjoint.cpp:1120-1132,62-64,2742; src/LaMEMLib.cpp
- Feedback requested from: Boris Kaus
- Status: open
- Boris: fix

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
- Introduced by (git blame at HEAD):
  - Garrett Ito — 53243845 DBDikeCreate:  creation of DA_CELL_1D and DA_CELL_2D_tave done ONCE he (2023-03-22); 2847e43e Dike.cpp Smooth_Sxx_eff now does time averaging over istep_nave steps  (2023-01-13); 8ac04ecf dike.cpp Communication arrays ycoords, ybounds, and sxx working for ne (2023-04-28) — DA_CELL_1D/2D_tave created with default partition, sized from first dike
  - Garrett Ito — c2f25559 Updates since last sync with bitbucket/lamem kaus dike.cpp:  Gaussian  (2023-08-06); e58ddd84 dike.cpp 2D elliptical Gaussian weighted filter on gsxx_eff_ave workin (2023-04-30) — y-neighbour exchange and istep_nave check without up-front limit checks
- At HEAD: still present, src/dike.cpp:108-120,959-1001,1314,1397
- Feedback requested from: Garrett Ito
- Status: open

## Notes

### N-001 [note, confirmed] `dikeRHS` is not set when dikes are off
- Tags: note, confirmed, logic, quick
- Where: src/JacRes.cpp:917, passed at :1079
- What: it's only read inside cellConstEq when `actDike` is on, so it's safe today.
  Compilers warn about it, and any change to the `if` order would read garbage.
- Fix idea: initialize it to 0 at the declaration.
- Introduced by (git blame at HEAD):
  - Jana Schierjott — 466d73a9 moved the dike RHS function to constEq.cpp in the end. the rest is the same exce (2021-04-28) — declared dikeRHS; assigned only inside actDike branch
- At HEAD: still present, src/JacRes.cpp:917 (decl), 1008-1011, 1079
- Feedback requested from: Garrett Ito
- Status: open

### N-002 [note, confirmed] Material data changed inside compute loops; Controls copied per cell
- Tags: note, confirmed, cleanup, quick
- Where: src/JacResAux.cpp:389-390 (`mat->rp` clamped), src/JacResTemp.cpp:107-110
  (`M->nu_k = 1`), src/JacResTemp.cpp:85 (`Controls ctrl = jr->ctrl` copied on every call)
- Why it matters: input checks belong in setup. The struct copy runs for every cell on every
  residual evaluation, which costs time for no reason.
- Introduced by (git blame at HEAD):
  - Jana Schierjott — bf918d9a changes for the temperature-dep conductivity, renamed the input to nu_k and move (2021-07-08); 6c30669d new way for the conducitvity seems to work except for the conductivity is multip (2021-06-30) — ctrl=jr->ctrl copy per cell; M->nu_k=1.0 set in compute loop; declared Controls ctrl local in JacResGetTempParam
  - Anton Popov — 94177af5 . (2017-02-16) — mat->rp clamped inside pore-pressure loop
- At HEAD: still present, src/JacResTemp.cpp:67,85,107-110; src/JacResAux.cpp:389-390
- Feedback requested from: Anton Popov, Garrett Ito
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
- Introduced by (git blame at HEAD):
  - Anton Popov — 49845492 Stabilizations and regularizations. (2020-05-10); d3ad1f41 Viscoplastic regularization updates and fixes (2024-02-21) — added eta_st stabilization in parallel to all phases; eta_min used as default eta_st (JacRes.cpp:150)
  - Boris Kaus — 4e46b0a7 fix to the previous commit (2020-05-20) — default eta_st = eta_min in phase.cpp global-override routine
- At HEAD: still present, src/JacRes.cpp:92,150; src/phase.cpp:1499; src/constEq.cpp:375,795-801,891-894
- Feedback requested from: Anton Popov, Boris Kaus
- Status: open

### N-004 [note, confirmed] Creep input checks never fire; `Bd`/`Bn` silently overwritten
- Tags: note, confirmed, input, quick
- Where: src/phase.cpp:526-548
- What: the diffusion check lists all four eta/Bd combinations, so it can never fail. If both
  `eta` and `Bd` are given, `Bd` is replaced by `1/(2 eta)`. The dislocation check allows
  `eta0+e0+n` together with `Bn`, and `Bn` is then overwritten. The error messages claim
  the parameters must be unique.
- Fix idea: reject (eta && Bd) and (eta0 && Bn).
- Introduced by (git blame at HEAD):
  - Anton Popov — 5b1d45fe Automated code formatting using astyle (2026-07-06); 64d9ea3b . (2017-02-16) — check collapsed to always-true (eta&&Bd allowed) in 'formatting' commit; Bd/Bn silently overwritten from eta/eta0,e0
  - Panupong Kongpet — c8ec3931 Update LaMEMLib.cpp-remove PetscMemzero, but add conditions in phase.cpp for uni (2026-02-26) — added restart-compare alternative to uniqueness check
- At HEAD: still present, src/phase.cpp:514-523 (diffusion), 527-536 (dislocation). Note: Always-true check came from 5b1d45fe (nominal astyle commit, Anton) which replaced Panupong's tolerance comparison with bare (eta && Bd); overwrite logic is original.
- Feedback requested from: Anton Popov
- Status: open (moved from A to C: `eta && Bd` accepted on purpose for restarts)

### N-005 [note, confirmed] meltParam.cpp is dead code
- Tags: note, confirmed, cleanup, quick
- Where: src/meltParam.cpp (whole file)
- What: no function in it is called. The only reference is a commented-out line in
  constEq.cpp:144-147.
- Fix idea: remove it or wire it in, but don't leave it half-done.
- Introduced by (git blame at HEAD):
  - Jianfeng Yang ("Jeff") — 901668f1 add parameterized hydrous MANTLE melting by Katz et al. 2003, G3. a minor modifi (2020-09-23) — created Katz melt parametrization file (later renamed meltParam)
  - Andrea Piccolo — fb9f55d0 Merge branch 'Simplified_Phase_Transition' of https://bitbucket.org/bkaus/lamem  (2020-09-25) — commented-out Compute_Melt_Fraction call in constEq.cpp
- At HEAD: still present, src/meltParam.cpp, src/meltParam.h (whole); src/constEq.cpp:144-147 (commented call). Note: Jianfeng Yang appears as 'Jianfeng Yang <yangjf@mail.iggcas.ac.cn>' and 'Jeff'. File renamed by Jana d6ab40b8, Anton 410bc5de.
- Feedback requested from: Andrea Piccolo
- Status: open

### N-006 [note, confirmed] PrintMatProp prints wrong values under some labels
- Tags: note, confirmed, logic, quick
- Where: src/phase.cpp:1471, 1479
- What: it prints `rho_c` twice (`rho_n` is missing). In the plasticity line, `eta_st` is
  printed under the "ch" label and `eta_vp`/`eta_st` are swapped. Only used for debugging.
- Introduced by (git blame at HEAD):
  - Boris Kaus — 98381d8e added a routine with which we can print out all material properties of a phase ( (2020-05-20); 9b87e7ac initial version of viscoplastic regularisation (2024-02-12) — added PrintMatProp, rho_c printed twice; inserted eta_vp in format but args left eta_st,eta_vp (swapped)
- At HEAD: still present, src/phase.cpp:1458 (rho_c twice), 1467 (plasticity print)
- Feedback requested from: Boris Kaus
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
- Introduced by (git blame at HEAD):
  - Yun Luo — e4778c61 New coupling module based on LaMEM 3.0.1 (2026-08-02) — FastScape module: <=2 check, kd units/array comments, time_fs; dup block surf.cpp:243
  - Wenrong Cao — 749ae824 Add slope-dependent erosion to the free surface (v3.1.0) (2026-08-13) — unit_m_yr comment wrong in slope-erosion code
  - Tobias Baumann — 283d15db added another sedimentation process (2019-03-18) — marginE read twice in sedimentation model 2
- At HEAD: still present, src/surf.cpp:111-112,243-248,72-73,1275; src/fastscape.cpp:92,124-130,1528; src/fastscape.h:163. Note: Wenrong Cao email caowenrong@gmail.com; Tobias Baumann appears as 'tosaba'. FastScape author raw name 'Yun'.
- Feedback requested from: Wenrong Cao, Yun Luo, Anton Popov
- Status: fixed (group C): all items; FastScape folder time now (t+dt)*s

### N-008 [note, confirmed] Average topography is a plain node average
- Tags: note, confirmed, physics, quick
- Where: src/surf.cpp:288-292
- What: `avg_topo` is the mean over nodes, not weighted by area. On non-uniform grids,
  infinitely fast erosion (model 1) doesn't keep the volume.
- Introduced by (git blame at HEAD):
  - Anton Popov — 5a7234af salt domes setup + depth-dependent density (2015-05-27) — VecSum(gtopo)/(Nx*Ny*Nz) plain node average introduced
- At HEAD: still present, src/surf.cpp:288-292 (FreeSurfGetAvgTopo)
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Anton Popov — 896e6bfd ISA orientation computation algorithm after Kaminsky et al. (2004), GJI. The alg (2015-06-19); f053c7bf ISA debugged (2015-11-18) — ISA/eigen algorithm after Kaminsky 2004; ISA debugged
- At HEAD: still present, src/Tensor.cpp:107-310 (Tensor2RNEigen, getISA), 488-603 (Tensor2RSSpectral)
- Feedback requested from: Anton Popov
- Status: fixed (group C): removed

### N-010 [note, confirmed] Non-unit bias is rejected, so all bias code is dead
- Tags: note, confirmed, cleanup, quick
- Where: src/fdstag.cpp:80-83
- What: MeshSeg1DReadParam stops with "Non-unit bias ratios are deprecated". So the biased
  branches in MeshSeg1DGenCoord (fdstag.cpp:153-171), FSLoadNonUniformGrid and
  FastScapeCreateGlobalGrid (fastscape.cpp:578-592, 921-937) never run. `bias_x` is still
  read as an input.
- Fix idea: remove the dead branches, or drop the check if bias should come back.
- Introduced by (git blame at HEAD):
  - Anton Popov — 3ce9bc02 Removed outdated files, deprecated bias ratios (2024-10-23); 88b60a16 code integration & functionality restore update version 0.0 (2014-09-09) — deprecated non-unit bias (error), leaving biased branches dead; original bias mesh code
  - Yun Luo — e4778c61 New coupling module based on LaMEM 3.0.1 (2026-08-02) — FastScape grid copies biased branches, unreachable
- At HEAD: still present, src/fdstag.cpp:80-83; dead branches fdstag.cpp:144-171, src/fastscape.cpp:564-631,907
- Feedback requested from: Anton Popov, Yun Luo
- Status: fixed (group C): bias removed fully; `bias_*` keys no longer read

### N-011 [note, confirmed] Velocity cylinder/box edge cases
- Tags: note, confirmed, input, quick
- Where: src/bc.cpp:2405-2411, 2234, 2369
- What: with `vmag` and base == cap the axis length is 0, giving NaN velocities; radius 0
  divides by zero. No checks. Velocity boxes and cylinders are skipped during the initial
  guess solve, but Bezier blocks and boundary inflow are not.
- Introduced by (git blame at HEAD):
  - Arne Spang — df705292 Internal velcocity cylinders that are suitable to enforce flow along a pipe. Als (2022-07-11) — velocity cylinder: vmag*ax/a with no zero-length/radius check
  - Andrea Piccolo — 1bc3ab4d Testing the new boundary conditions, still working on it (2020-10-04) — initGuess skip for boxes/cylinders only
- At HEAD: still present, src/bc.cpp:2399-2411 (BCApplyVelCylinder vmag/a), 2234, 2369 (initGuess skip). Note: Arne Spang raw name 'Arne_Spang'.
- Feedback requested from: Andrea Piccolo, Anton Popov
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
- Introduced by (git blame at HEAD):
  - Andrea Piccolo — 199277db Fixing a bit of stuff. The function is there, but, it is not converging properly (2020-12-19); 920a2688 Bug fix: I need to clean the override phase function and introducing better cont (2020-11-12) — Influx_type message; plume code; hard-coded k_thermal in BCOverridePhase
  - Boris Kaus — b9d87321 * added a gaussian-type inflow/outflow velocity at the bottom, which is more smo (2020-10-20) — bc.h comment Gaussian=0/Poiseuille=1 wrong vs code
  - Panupong Kongpet — 7d0eb7a1 Add velin_net_num_periods, velin_net_time_delims, and bvel_velin_net for the net (2026-01-25) — BCGetVelins debug print with (Myr); velin_net duplication
- At HEAD: still present, src/bc.h:263; src/bc.cpp:782,827,1683,2961 (tSPCList partly used in JacResTemp.cpp:415). Note: Broad clean-up note; attributed to dominant authors. tSPCList sub-item is stale (used in JacResTemp.cpp).
- Feedback requested from: Andrea Piccolo, Boris Kaus
- Status: open

### N-013 [note, possible] AVD point with no cells gives NaN coordinates
- Tags: note, possible, logic, quick
- Where: AVD.cpp:1123-1125 (commit 6124aa7d)
- What: division by `hclaim` with no guard; `hclaim = 0` only if a point owns no AVD cells
  (unlikely).
- Introduced by (git blame at HEAD):
  - Adina Püsök — 766f54c3 1. Marker control (solves the "Empty Control Volume" error) * activate with '-us (2015-03-30) — AVD marker control: xc/hclaim with no zero guard
  - Anton Popov — 761a531d Let's be more accurate with the type conversion. (2015-07-27) — added (PetscScalar) cast on the division
- At HEAD: still present, src/AVD.cpp:1123-1125 (also 546-548). Note: Adina Puesoek raw name 'Adina Püsök'. Division by hclaim original from 766f54c3.
- Feedback requested from: Anton Popov
- Status: fixed (group A)

### N-014 [note, confirmed] Docs say Clapeyron parameters overwrite a predefined law
- Tags: note, confirmed, input, quick
- Where: info/options/input_file.dat:743-749; src/phase_transition.cpp:562-580
  (commit 6124aa7d)
- What: the example sets `Name_Clapeyron = Mantle_Transition_660km` together with
  `clapeyron_slope`, `P0_clapeyron`, `T0_clapeyron` marked "overwrite". The code reads
  these only for `user-defined`, so they are ignored.
- Fix idea: read them as optional overrides after `SetClapeyron_Eq`, or fix the docs.
- Introduced by (git blame at HEAD):
  - Anton Popov — 326a1b51 Finally clear logic in Set_Clapeyron_Phase_Transition (2026-04-29) — made clapeyron_slope/P0/T0 read only for user-defined; overrides dropped
  - Arne Spang — bf933543 Added explanation how to overwrite parameters (2025-02-24) — docs say these parameters 'overwrite' predefined law
- At HEAD: still present, info/options/input_file.dat:743-746; src/phase_transition.cpp:560-580. Note: Arne's earlier code (bf933543) honored overrides; Anton's 2026-04-29 rewrite removed them, docs stayed.
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Tobias Baumann — a6a5e461 Added objective function (2015-06-15) — created objective function module (ObjFunct)
  - Georg Reuber — 625d00d2 added FD kernel calculation for eta and n (2020-06-06) — added FD kernel routines cellConstEqFD/edgeConstEqFD/devConstEqFD/setUpPhaseFD
- At HEAD: still present, src/objFunct.cpp (all); src/adjoint.cpp:3684-4150. Note: Tobias Baumann appears as 'tosaba'; Georg Reuber as 'greuber'.
- Feedback requested from: Anton Popov, Boris Kaus
- Status: fixed (group C): objFunct.cpp, FD kernels and unused fields removed

### N-016 [note, confirmed] `gw_level_type = surf` gives a flat water table
- Tags: note, confirmed, physics, medium
- Where: src/JacResAux.cpp:357
- What: it uses `surf->avg_topo`, the mean topography, not the local surface height.
- Why it matters: with real topography, pore pressure is too high under valleys and too low
  under hills. Users may expect the water table to follow the surface.
- Introduced by (git blame at HEAD):
  - Anton Popov — 94177af5 . (2017-02-16) — gwLevel = surf->avg_topo for gw_level_type surf
- At HEAD: still present, src/JacResAux.cpp:357
- Feedback requested from: Anton Popov
- Status: open

### N-017 [note, possible] Creep profile experiment types not checked
- Tags: note, possible, physics, medium
- Where: src/phase.cpp:888-1352
- What: each profile sets `_UniAxial_` or `_SimpleShear_`, and that changes the prefactor
  by up to 3^((n+1)/2)/2 vs 2^(n-1). For example, the HK2003 laws are marked simple shear,
  but most of their data are triaxial. Needs a check against each paper.
- Introduced by (git blame at HEAD):
  - Anton Popov — cb50cfeb dc&ps creep (started) (2019-06-24); 64d9ea3b . (2017-02-16) — assigned _SimpleShear_/_UniAxial_ types to creep profiles; original creep profile database
  - Boris Kaus — aef1985b added dry plagioclase rheology (diffusion & dislocation creep). (2020-10-15) — added plagioclase profiles with experiment types
- At HEAD: cannot be judged by blame, src/phase.cpp:840-1352 (SetDiffProfile/SetDislProfile). Note: Possible physics issue, not verifiable by git; attribution by who set profile types.
- Feedback requested from: Anton Popov, Boris Kaus
- Status: not an issue (group C, Anton)

### N-018 [note, confirmed] Polygon helpers use stack arrays and can divide by zero
- Tags: note, confirmed, memory, medium
- Where: src/tools.cpp:608-609 (interpStretch), 658-662 (findCenterMass); marker.cpp:1676-1677
- What: arrays sized by the number of polygon nodes/slices live on the stack (variable-length
  arrays, not standard C++). Large polygons can overflow the stack. findCenterMass divides by
  the area with no check, so a degenerate polygon gives inf/NaN.
- Introduced by (git blame at HEAD):
  - Arne Spang — a88df698 added geometry manipulation (2020-07-08); ec141cb8 The old version worked on pluton but not on mogon (2020-08-26); 82be9435 Extended ability to read geomIO polygons (2022-05-05) — findCenterMass VLAs, no zero-area check; interpStretch VLAs; SxAll/SyAll VLAs in marker.cpp
- At HEAD: still present, src/tools.cpp:608-609,658-662,704-705; src/marker.cpp:1676-1677
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Boris Kaus — e173ff0e improved the manner in which we switch from Newton back to Picard and viceversa (2015-03-20) — removed refRes=nrm at switch; divergence test vs first-iteration residual
  - Anton Popov — ab24f955 Jacobian type control (sory for delay, I wasn't actually doing that) (2015-03-12); b7f5b864 default solver settings (update) (2025-08-22) — original Newton/Picard switch and 'MMFD' label; maxItNwt/rtolNwt defaults
  - Arne Spang — 00b99981 fixed solver option -snes_NewtonSwitchToPicard_it (2020-04-02) — itNwt counter: Newton-iteration cap then Picard step
- At HEAD: still present, src/nlsolve.cpp:209-232 (switch logic), 247 (MMFD label)
- Feedback requested from: Anton Popov, Boris Kaus
- Status: open

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
- Introduced by (git blame at HEAD):
  - Anton Popov — ccb11b49 restructuring (explanation follows) (2014-10-24); ca1e35aa Status quo but with new data structures (2025-04-23); cdf449a2 ... (2025-04-23) — original pgamma / inv_eta Schur code; jp_type user/bf handling; pgamma checks in lsolve
- At HEAD: still present, src/lsolve.cpp:48-90; src/matrix.cpp:316,1080,1131,1147
- Feedback requested from: Anton Popov
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
- Introduced by (git blame at HEAD):
  - Anton Popov — e5372907 SNES part III (2014-09-24); 71073db8 Added clean assembly functions need by PCMG (2025-04-25); 0b793ef0 update ... (2025-04-20) — SNES/KSP option ordering, PCMAT; unused multigrid test helpers; MGLevelCreate 3-arg declaration
  - Boris Kaus — 199c46dd my last commit didn't work so well with many of the tests (apologies!). I now ro (2021-02-03) — MatMFFDSetType(WP) overrides -fd_mat_mffd_type
- At HEAD: still present, src/nlsolve.cpp:79-93,286; src/multigrid.cpp:283-284,1603-1802; src/multigrid.h:88; src/matData.cpp:159; src/lsolve.cpp:82-90; src/matAux.cpp:65,124-127. Note: Broad clean-up note; nearly all Anton.
- Feedback requested from: Anton Popov, Boris Kaus
- Status: partly fixed (group C): dead multigrid test functions removed, MGLevelCreate declaration fixed; rest open

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
- Introduced by (git blame at HEAD):
  - Adina Püsök — 766f54c3 1. Marker control (solves the "Empty Control Volume" error) * activate with '-us (2015-03-30); f848a86c CVI: *stand-alone advection schemes and velocity interpolation routines *mainly  (2015-07-20) — AVD dead injection functions, buffer=1; CVI 'RK 4th order' comment
  - Anton Popov — c630c2e4 subgrid finished (2019-04-16); c991f16a Get rid of PetscMPIInt in LaMEM stuctures (2026-06-13) — subgrid test block; nbyte cast to PetscMPIInt
  - Boris Kaus — 8dd159de added possibility to perturb the geometrical object 'layers' with an additional  (2020-05-22) — marker.cpp layer comment says sin, code uses cos
- At HEAD: still present, src/AVD.cpp:98,403-665; src/advect.h:71; src/subgrid.cpp:740-786; src/cvi.cpp:13,157; src/marker.cpp:974,1071-1082; advect.cpp:1137; cvi.cpp:683,696. Note: Broad note. Others: Jianfeng Yang (_DISP_ adfab528), Arne Spang (ellipsoid parser 28afee1c lacks fb->ID).
- Feedback requested from: Anton Popov, Boris Kaus
- Status: partly fixed (group C): dead AVD functions, `_DISP_`, subgrid test block removed; rest open

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
- Introduced by (git blame at HEAD):
  - Garrett Ito — 603d8147 phase_transition: fixed bug in Check_NotInAirBox...P->X[2] was being used to int (2022-10-15); eaecfa33 Updates to dike.cpp, t26_dike, and lamem_input.dat since 8/23 pull request to Ma (2025-10-21); 1592aab1 dike.cpp:  Using module function in dynamic diking routines, cleaned up some lin (2023-08-09) — istep+1==0 dead test; debugging condition left in Locate_Dike_Zones; magma_presence=1.0 //testing
  - Andrea Piccolo — d7aff650 Bug Fix: We need to communicate to LaMEM the case in which no phase change are r (2022-03-02) — PhaseInside[0]=-1 set after 'Don't set phase' message
  - Jana Schierjott — 6b60a55b added new parameters for the location of Mc (2021-11-17) — unused x_Mc,z_Mc in dike.h
- At HEAD: still present, src/dike.cpp:508,522-523,748,1437; src/dike.h:44-46; src/phase_transition.cpp:130,333,451,975,1478-1515; phase_transition.h:49. Note: Broad note; Garrett Ito raw 'Garrett Ito'/'Garrett'. Piccolo raw 'andrea <andrea@Piccolo-W10.bgi.uni-bayreuth.de>'.
- Feedback requested from: Garrett Ito, Andrea Piccolo
- Status: open
- Boris: Garrett should double check. Pressure is always in Pa, also with GEO units, so that is fine. Most others seem good to be improved by GI should check.

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
- Introduced by (git blame at HEAD):
  - Anton Popov — 722900e1 debug (2017-02-28); 0627e8ff dissolved user context and added new input framework (2016-11-19) — strtol/strtod unchecked parsing; new input framework (FBLoad, rank-0 errors)
- At HEAD: still present, src/parsing.cpp:370,431 (strtol/strtod), 507-570, FBLoad ~56-87
- Feedback requested from: Anton Popov
- Status: fixed (group B): approved parts; parser relaxed
- Boris: many of these make sense and would make the code more robust

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
- Introduced by (git blame at HEAD):
  - Anton Popov — 229a37ff Automatic tolerances polished (2026-04-23); b7f5b864 default solver settings (update) (2025-08-22); a912ef63 scaling, time step settings, main framework, input file (2016-11-20) — its_ steady-state temp solver shares ts_ksp_ref_norm; %g option storing 6 digits; scaling labels / typos
  - Boris Kaus — 35792000 This adds the possibility to not write any output  at all by specifying -nstep_o (2023-06-02); c6cb8df4 add documentation to /docs directory in preparation to use Documenter to create  (2023-05-11) — nstep_out<=0 disables output; docs say nstep_out sets number of steps
  - Arne Spang — 1c27c92e New options that allow planning the timesteps for the whole simulation (2022-05-04) — fix_dt schedule applied one step late
- At HEAD: still present, src/tssolve.cpp:216,234-296; src/options.cpp:906-907; src/scaling.cpp:118,183,241; src/LaMEMLib.cpp:986-1026; doc/src/man/GettingStarted.md:49-60
- Feedback requested from: Anton Popov, Boris Kaus
- Status: partly fixed (group B): labels/typos fixed; fix_dt, dt_out, time_end trimming left for C
- Boris: the obvious ones can be fixed; nstep<=0 should indeed switch off output

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
- Introduced by (git blame at HEAD):
  - Anton Popov — b1d294ea Phase aggregate output (done) (2019-03-19); ade15fb9 * pvd file is updated optionally with every output, controlled by option:    -ou (2014-11-04); 4d4bde66 [1] Geometric primitives are applied in the order of appearance in the input fil (2017-06-28) — phase aggregate output without negative-ID check; stub output functions; use_bound handling for total_pressure
  - Jianfeng Yang ("Jeff") — 6c7c0b5e Output total dissipation rate (shear heating). effective shear heating is only c (2016-04-09) — plast_dissip outputs total shear heating Hr
  - Marcel Thielmann — 9fa51713 it now outputs the energy residual if the option -out_energ_res 1 is specified o (2015-12-02) — energ_res output of stale jr->ge
- At HEAD: still present, src/outFunct.cpp:131,429-488,617-640,672; src/paraViewOutBin.cpp:286; paraViewOutAVD.cpp:365,909; paraViewOutPassiveTracers.cpp:319. Note: Broad note; Jianfeng Yang raw 'Jeff <jianfengyang1989@gmail.com>'.
- Feedback requested from: Anton Popov
- Status: partly fixed (group B): energ_res, pShift, AVD, outputs fixed; fwrite per value and total_pressure edges left
- Boris: plast_dissip is correct; energ_res to be fixed; vol_rate, vorticity, ang_vel_mag to get an option to write; Passive tracer pressure should get p_shift; out_avd_ref <= 0 should be rejecetd

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
- Introduced by (git blame at HEAD):
  - Boris Kaus — 3751f5eb more rewriting (2020-05-13); 85c78906 WIP: added routine to change the LaMEM material database (2020-05-15); fec8c737 set boundary conditions for the adjoint system (2020-06-10) — optional-but-required Gradient/Scale params; stack Scaling pointer in modified DB; Adjoint_ApplyBCs sets dF to BC values
  - Georg Reuber — cf0ea953 Added the adjoint gradient computation and the TAO framework to LaMEM as wrapper (2018-01-11); fdbed43e Clean up of LaMEM.cpp. (2018-01-15) — AdjointPointInPro and observation points; !IOparam->Ap check
- At HEAD: still present, src/adjoint.cpp:449-458,804,903,951,1609,1701,1768-1776,2184-2199,2879,4155-4187. Note: Georg Reuber raw 'greuber/Georg Reuber <georg.reuber@gmx.de>'.
- Feedback requested from: Boris Kaus
- Status: open
- Boris: try fixing and see if tests still work

### N-028 [note, confirmed] FSSA only stabilizes the diagonal terms
- Tags: note, confirmed, physics, design
- Where: src/JacRes.cpp:1086-1103
- What: each momentum equation gets only its own velocity times its own density gradient
  (for example vz*drho/dz*gz). The full term is (v.grad rho)*g, which would include
  vx*drho/dx*gz.
- Why it matters: lateral density gradients (sloped surfaces, tilted interfaces) are not
  stabilized. The code comment says "lumped", but the user docs don't mention the limit.
- Introduced by (git blame at HEAD):
  - Anton Popov — 40af2e14 Density gradients stabilization terms are added to the residual and precondition (2015-02-25) — FSSA density-gradient stabilization, diagonal terms only
  - Boris Kaus — 199c46dd my last commit didn't work so well with many of the tests (apologies!). I now ro (2021-02-03) — restored diagonal-only form after trying full v-sum
- At HEAD: still present, src/JacRes.cpp:1086-1103. Note: Boris tried (vx+vy+vz)*t then rolled back to own-component form in 199c46dd.
- Feedback requested from: Anton Popov, Boris Kaus
- Status: won't fix (group C): Boris: leave
- Boris: leave as-is, as this is the key contribution. I'm pretty sure I had an option to do the full terms @ some stage

### N-029 [note, confirmed] Rock/air phase change keeps old stress, APS and T
- Tags: note, confirmed, physics, design
- Where: subgrid.cpp:528-606 (commit 6124aa7d)
- What: `ADVMarkCrossFreeSurf` changes only the phase; history is mixed into cut cells for one
  step. There is already a WARNING comment.
- Introduced by (git blame at HEAD):
  - Anton Popov — d3a5c7ee marker-surface-interaction major overhaul (2019-04-02); 60dc961e marker-free-surface-interaction update (2019-03-26) — marker-surface interaction overhaul; phase-only change, WARNING comments; initial marker/free-surface rock-air logic
- At HEAD: still present, src/subgrid.cpp:537-615 (ADVMarkCrossFreeSurf)
- Feedback requested from: Anton Popov
- Status: won't fix (group C): Boris: leave
- Boris: don't think its easy to figure out what that old stress should be if we go from air->rock (other way is easier), so I would leave as-is
