# dylib_plugins

Julia sources for LaMEM's user-defined dylib plugin (`src/dylib_plugins.h`).
A rule file `include`s `LaMEMPlugin.jl` and defines
`rule(m::MarkerView) -> (phase, T_dimensional)`, wrapped by
`lamem_pt_wrapper` into the `lamem_phase_transition`/
`lamem_phase_transition_abi_version` `@ccallable`s that `-dylib_plugin`
loads. `ptlib_constant.jl`/`ptlib_box.jl` are example rules, one per file.

Build: `julia --project=@juliac build_plugin.jl <rule.jl> [output_dir]`
(`@juliac`: an environment with `JuliaC` installed). Produces
`<output_dir>/build_<name>/lib/libptlib_<name>.{dylib,so}`.

.dat keywords: `dylib_plugin = <path>` loads the library once at startup;
`phase_transitions = builtin | dylib` selects whether
`lamem_phase_transition` also runs each step.

Only one plugin library loads per LaMEM process. It embeds its own Julia
runtime and only works as a standalone LaMEM run - it cannot be `ccall`ed
from inside an already-running Julia process. CI builds/tests with
Julia 1.13 + JuliaC.jl. Windows is not yet supported (dylib_plugin fails
with a clear error there).
