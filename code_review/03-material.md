# Task 3: Material laws

Commit: 4835b6ba   Started: 2026-10-16   Status: done (first pass)

## Coverage
- [x] constEq.cpp / constEq.h: all functions
- [x] phase.cpp: softening, phase reading, checks, scaling, creep profiles, unit corrections
- [x] meltParam.cpp: only checked for callers. It's dead code (N-005), so the formulas weren't reviewed.
- matData.cpp moved to task 4 (it's matrix setup, not material data).
- Phase diagram file loader (marker.cpp:1740-1840) read to check S-031.

## Notes and open questions
- Mixed cells average phase viscosities arithmetically (constEq.cpp:525). That's a design
  choice, but a strong phase dominates a mixed cell. Worth knowing when comparing to codes
  that use harmonic or geometric averaging.
- `p_lim_plast` is switched off after the first iteration in nlsolve.cpp:254, so it does
  work as documented.
- Simple-shear correction 2^(n-1) and uniaxial 3^((n+1)/2)/2 match Gerya (2010), assuming
  the lab law is written with differential stress. Whether each profile's experiment type
  (for example HK2003 = simple shear) is right needs a check against the papers (N-017).
- Fracture-energy regularization (`Lm`) has no documentation anywhere (S-002).

## Issues
Moved to findings.md as S-031, M-003, S-002, S-003, M-048, N-003, N-004, M-004, M-049, N-005, N-006, N-017.
