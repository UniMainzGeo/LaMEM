# Task 10: Input, time stepping, driver, passive tracers

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] parsing.h / parsing.cpp: all functions
- [x] options.h / options.cpp: all functions
- [x] scaling.h / scaling.cpp: all functions
- [x] tssolve.h / tssolve.cpp: all functions
- [x] LaMEMLib.h / LaMEMLib.cpp: all functions (create, restart save/load, time loop, temperature init)
- [x] LaMEM.h / LaMEM.cpp: main
- [x] passive_tracer.h / passive_tracer.cpp: all functions
- Checked call sites: advect.cpp (ADVExchange, ADVRemap, ADVMapMarkToCells, ADVReAllocStorage),
  tools.cpp (Dir* helpers), fdstag.cpp (FDSTAGGetLevelsLocalGridSize), phase_transition.cpp (MovingBox)
- Not run: all findings are from reading the code.

## Notes and open questions
- Start-up order: FBLoad -> setSolverOptions (clears PETSc options, sets defaults, reads
  <PetscOptions> block, then re-inserts the saved command-line options at the end) -> LaMEMLibMain.
  Everything read while the options are cleared ignores the command line (S-039).
- PETSc gives priority to the option inserted last, so the re-insert step does make command-line
  options win for the solver settings themselves.
- Block keys on the command line are `-key[ID]` with `fb->ID`. For the SolverOptions block `ID` is
  whatever was left from before; it doesn't matter because options are cleared there anyway.
- No block keyword is a substring of another (checked `<BoxStart>` vs `<VelBoxStart>` etc.).
- Whole LaMEMLib struct (TSSol, Scaling, P_Tr settings) is fwrite'd to the restart file, so the
  input file has no effect on those after restart (except DBMat and dikes, see M-034).
- Scaling: GEO conversions checked (Myr, km, cm/yr, MPa, mW/m^2, Tshift = 273.15). All right.
- TSSolIsOutput has a side effect (updates time_out) but it's called once per step only.
- Output with `nstep_out <= 0` is off completely (also the initial step and dt_out).
- Passive tracers: every rank holds all tracers; ranks that don't own a tracer write -DBL_MAX and
  an Allreduce(MAX) merges them. A tracer that is in no local box (outside the domain, or exactly
  on the global upper bound) gets -DBL_MAX for good.
- Tracer advection is forward Euler with the end-of-step velocity, done after BCStretchGrid
  (grid already stretched), while markers can use RK schemes. Tracers and markers drift apart.

## Issues
Moved to findings.md as S-024, S-038, S-039, M-036, M-066, M-037, M-038, M-067, M-068, N-024, N-025.
