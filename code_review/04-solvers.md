# Task 4: Solvers

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] nlsolve.cpp / nlsolve.h: all functions (setup, Picard/Newton switch, MFFD, line-search
      pre-check, convergence test, failure report, automatic tolerances)
- [x] lsolve.cpp / lsolve.h: all functions (option parsing, coupled MG, block factorization,
      wBFBT apply, user fieldsplit)
- [x] matData.cpp / matData.h: all functions (index vectors, parameter setup, coarsening of
      parameters and BCs, SPC lists)
- [x] matrix.cpp / matrix.h: all functions (preallocation, monolithic and block assembly,
      Picard products)
- [x] matAux.cpp: all functions (AIJ helpers, near null space, local stiffness, constraints,
      block scatter)
- [x] matBFBT.cpp: all functions
- [x] matFree.cpp / matFree.h: all functions (operator, diagonal, restriction, prolongation)
- [x] multigrid.cpp / multigrid.h: all functions; 2D transfer operators only skimmed
- Not checked: the sign and scale of every stencil entry against JacResFormResidual
  (spot-checked x-direction normal and xy shear terms, pressure gradient and FSSA; they match).

## Notes and open questions
- Assembled Picard rows for constrained DOFs are identity; matrix-free rows are zero (diagonal
  reports 1). Works because residuals at constrained DOFs are always zero, but the two
  operators differ. Added to N-020.
- With pgamma = 1 the monolithic preconditioner still has -1/eta on the pressure diagonal
  (a stabilization). Picard product subtracts it again via the M matrix. Consistent.
- Block upper/lower factorization and wBFBT (K = Apv*C*Avp, C = 1/sqrt(eta_face)) checked
  against the assembled blocks: consistent, signs right.
- Galerkin restriction is not a scaled transpose of prolongation (transverse injection vs
  bilinear), so coarse operators are not symmetric. Weights assume a uniform grid. Both are
  valid choices for a preconditioner; noted only.
- MFFD residual calls overwrite svCell/svEdge with the perturbed state. The temperature
  pre-check then uses this state. The perturbation is ~1e-8, so no practical effect.
- `MatDataCoarsen` doesn't copy `rescal`. Harmless: `rescal` is only allowed with Galerkin
  coarse levels, which never use the coarse MatData for the operator.
- Test coverage: no test uses `periodic`, and none uses `-gmg_mat_free_levels`. Only t27 uses
  `use_mat_free_jac`. M-055, M-017, M-056 live in these untested paths.

## Issues
Moved to findings.md as M-055, M-017, M-056, M-018, M-019, N-019, N-020, N-021.
