# CLAUDE.md

## Build

- Build LaMEM in both modes after every change: `make mode=deb surf=scape all` and
  `make mode=opt surf=scape all` in `src/` (`surf=scape` enables FastScape, needed by t37;
  run `clean_all` after changing it). The tests use both binaries (most use opt, a few use
  deb), and hand runs of `bin/opt/LaMEM` or `bin/deb/LaMEM` with a stale build give wrong
  conclusions.
- If linking fails with undefined symbols, the objects are stale (built against another PETSc):
  run `make mode=deb clean_all` (or `mode=opt`) in `src/`, then rebuild.

## Tests

- Tests live in `test/` and run through Julia (`make test`).
- Use test selectors and run only the tests that touch the changed code: `make test 30 31` or
  `make test 03-07`. The full suite takes more than 10 minutes; run it only when needed.
- `make update NN` rewrites all reference (`.expected`) files in those test folders. Keep them
  all; references of tests that passed before do the same job. Say in the commit when a change
  alters results.

## Code review

- The review of `src/` is in `code_review/` on the `code-review` branch (keep it off `master`).
  Start with `code_review/README.md`; each reviewer's fix plan is `code_review/fix_plan_*.md`.
- Check each finding in the code before fixing it; some fix ideas are wrong.
