# dylib_plugins

User documentation: `doc/src/man/JuliaPlugins.md` ("Julia phase-transition
plugins" in the LaMEM manual). This file is the short developer summary.

Julia sources for LaMEM's user-defined dylib plugin, ABI v3
(`src/dylib_plugins.h`, mirrored struct-for-struct in `LaMEMPlugin.jl`).
A rule file `include`s `LaMEMPlugin.jl`, defines
`rule(m::MarkerView) -> MarkerView` (return `m`, or `update(m; phase=..,
T=.., aps=.., ...)`; a `(phase, T)` tuple also works) and one
`lamem_phase_transition` `@ccallable` that calls
`lamem_pt_wrapper(rule, markers, cells, step, scaling)`. `LaMEMPlugin.jl`
itself exports the `lamem_plugin_abi_version`/`lamem_plugin_struct_sizes`
`@ccallable`s LaMEM checks at load time.

`MarkerView` carries every marker field (dimensional); writable are phase,
T, p, APS, ATS, the deviatoric stress and the displacement, read-only is the
position. The marker p is the pressure history (projected to the cells as
`pn`, enters the next solve via `-IKdt*(p - pn)`), so changing it only acts
on compressible phases. `m.cell` is a lazy, dimensional view of all
`SolVarCell` fields of the marker's cell plus the cell-centred J2
invariants; `phase_ratio(m.cell, ph)` gives the phase ratios.

Examples, one rule per file: `ptlib_constant.jl`/`ptlib_box.jl` reproduce
the built-in Constant/Box transitions bit-for-bit, `ptlib_demo_fields.jl`
reads cell data and writes APS and p.

Build: `julia --project=@juliac build_plugin.jl <rule.jl> [output_dir]`
(`@juliac`: an environment with `JuliaC` installed). Produces
`<output_dir>/build_<name>/lib/libptlib_<name>.{dylib,so}`.

.dat keywords: `dylib_plugin = <path>` loads the library once at startup;
`phase_transitions = builtin | dylib` selects whether
`lamem_phase_transition` also runs each step.

Only one plugin library loads per LaMEM process. It embeds its own Julia
runtime and only works as a standalone LaMEM run - it cannot be `ccall`ed
from inside an already-running Julia process. CI builds/tests with
Julia 1.13 + JuliaC.jl, including on Windows.
