# User-defined phase transitions with Julia plugins

LaMEM's built-in phase transitions (`Constant`, `Clapeyron`, `Box`, see the `<PhaseTransitionStart>` blocks in `input_file.dat`) cover the common cases. When you need something they cannot express — a rule that depends on several marker quantities at once, a reaction with its own kinetics, a lookup table, a time-dependent threshold — you can write the rule in Julia instead and have LaMEM call it.

The rule is compiled once, ahead of time, with [JuliaC.jl](https://github.com/JuliaLang/JuliaC.jl) into a self-contained shared library (a *plugin*). LaMEM loads that library at startup and calls it once per time step, on every MPI rank, with all the markers that rank owns. There is no Julia process involved at run time and no interpreter overhead: the plugin is native code.

!!! note
    A plugin embeds its own Julia runtime. It works whenever LaMEM runs as its own process — including when it is launched from Julia through `LaMEM.jl`'s `run_lamem`, which starts LaMEM as a subprocess. It cannot be used from a LaMEM that is itself loaded as a library into an already-running Julia session. Only one plugin can be loaded per LaMEM run.

Plugins are supported on Linux, macOS and Windows, and are tested in CI on all three. You need Julia 1.13 or later and the `JuliaC` package.

## Writing a rule

A plugin is a single Julia file. It `include`s the helper module `scripts/dylib_plugins/LaMEMPlugin.jl` from the LaMEM repository, defines one function — the rule — and ends with a few lines of fixed boilerplate that export the rule to LaMEM.

The rule receives one marker at a time as a `MarkerView` and returns the marker's new phase and new temperature:

```julia
module PTLibMyRule

include(joinpath(@__DIR__, "LaMEMPlugin.jl"))   # adjust the path to where you keep it
using .LaMEMPlugin

const PHASE_BELOW = Cint(2)
const PHASE_ABOVE = Cint(3)

# Called once per marker. Return (new_phase, new_T); return (m.phase, m.T)
# to leave the marker untouched.
function my_rule(m::MarkerView)
    if m.phase == PHASE_BELOW || m.phase == PHASE_ABOVE
        new_phase = m.T >= 1200.0 ? PHASE_ABOVE : PHASE_BELOW
        return (new_phase, m.T)
    end
    return (m.phase, m.T)
end

# --- boilerplate, identical in every plugin ---------------------------------
Base.@ccallable function lamem_phase_transition(
        n::Csize_t,
        x::Ptr{Cdouble}, y::Ptr{Cdouble}, z::Ptr{Cdouble},
        T::Ptr{Cdouble}, p::Ptr{Cdouble}, time::Cdouble,
        sxx::Ptr{Cdouble}, syy::Ptr{Cdouble}, szz::Ptr{Cdouble},
        sxy::Ptr{Cdouble}, sxz::Ptr{Cdouble}, syz::Ptr{Cdouble},
        j2_stress_cell::Ptr{Cdouble}, j2_strainrate_cell::Ptr{Cdouble},
        eta_cell::Ptr{Cdouble}, aps_cell::Ptr{Cdouble},
        phase_in::Ptr{Cint}, phase_out::Ptr{Cint},
        T_out::Ptr{Cdouble},
        scaling::Ptr{LaMEMPluginScaling})::Cint
    try
        return lamem_pt_wrapper(my_rule, n, x, y, z, T, p, time,
            sxx, syy, szz, sxy, sxz, syz,
            j2_stress_cell, j2_strainrate_cell, eta_cell, aps_cell,
            phase_in, phase_out, T_out, scaling)
    catch
        return Cint(-1)
    end
end

Base.@ccallable function lamem_plugin_abi_version()::Cint
    return LaMEMPlugin.ABI_VERSION
end

end # module
```

`lamem_pt_wrapper` does all the bookkeeping: it converts LaMEM's internal non-dimensional arrays to dimensional values, builds a `MarkerView` for each marker, calls your rule, and writes the results back. Your rule only ever deals with physical quantities.

### What the rule sees

All `MarkerView` fields are **dimensional**, in the same units as your `.dat` file (for `units = geo`: km, Myr, °C, MPa, Pa·s).

| Field | Meaning |
|:--|:--|
| `x`, `y`, `z` | marker position |
| `T` | temperature |
| `p` | pressure (the same pressure the rheology sees) |
| `time`, `dt`, `step` | current time, time step size, step number |
| `sxx`, `syy`, `szz`, `sxy`, `sxz`, `syz` | deviatoric stress components |
| `j2_stress_cell`, `j2_strainrate_cell` | second invariants of stress and strain rate, taken from the marker's cell |
| `eta_cell` | effective viscosity of the marker's cell |
| `aps_cell` | accumulated plastic strain of the marker's cell |
| `phase` | current phase ID |
| `T_internal` | temperature in LaMEM's internal units, see below |

Cell-based quantities (`*_cell`) are the values of the finite-difference cell the marker sits in; they are not interpolated to the marker position.

### What the rule may change

Only the **phase** and the **temperature**. Everything else is read-only. A plugin cannot, for instance, reset the accumulated plastic strain the way the built-in `ResetParam = APS` option does.

Return the input temperature `m.T` unchanged whenever your rule does not modify it. The wrapper only converts a temperature back to internal units if it actually changed, so an unchanged `T` is passed through bit-for-bit; recomputing an "unchanged" value could perturb it by a rounding error.

### Constraints inside the rule

The plugin is compiled with `--trim=safe`, which keeps the library small but restricts what the rule may do:

* **No printing or other I/O** inside the rule. It is not supported in trimmed code and the plugin will fail to build.
* **Keep it type-stable.** Avoid dynamic dispatch, containers of abstract type, and global non-`const` state. Precompute constants with `const` at module level, as in the example.
* **Exceptions abort the run.** The boilerplate catches them and reports failure to LaMEM, which then stops with an error. Do not rely on exceptions for control flow.

### Reproducing a built-in transition exactly

LaMEM compares temperatures in its internal non-dimensional units. If you need a rule that reproduces a built-in transition bit-for-bit (for example to validate a plugin against a built-in run), compare `m.T_internal` against a threshold you non-dimensionalise with the same formula LaMEM uses, which `LaMEMPlugin` provides as `nondimensionalize_T`. `scripts/dylib_plugins/ptlib_constant.jl` does exactly this and documents the reasoning.

## Building the plugin

Install `JuliaC` into a named environment once:
```
julia --project=@juliac -e 'using Pkg; Pkg.add("JuliaC")'
```

Then build with the helper script from the LaMEM repository:
```
julia --project=@juliac scripts/dylib_plugins/build_plugin.jl ptlib_myrule.jl [output_dir]
```

This compiles the rule together with the Julia runtime it needs into a *bundle* directory. The plugin name is the file name without the extension and without a leading `ptlib_`, so `ptlib_myrule.jl` becomes:

| Platform | Library |
|:--|:--|
| Linux | `<output_dir>/build_myrule/lib/libptlib_myrule.so` |
| macOS | `<output_dir>/build_myrule/lib/libptlib_myrule.dylib` |
| Windows | `<output_dir>/build_myrule/bin/libptlib_myrule.dll` |

`output_dir` defaults to the directory containing the rule file. Building takes a few minutes, as the Julia runtime is compiled into the bundle.

!!! warning
    The library depends on the other files in its bundle directory (`libjulia` and friends). Move or copy the whole `build_<name>` directory, never the library file alone, or LaMEM will not be able to load it.

The two example rules in `scripts/dylib_plugins/` build the same way; `ptlib_constant.jl` and `ptlib_box.jl` reproduce the built-in `Constant` and `Box` transitions.

## Using the plugin in a model

Two keywords in the `.dat` file activate a plugin:
```
    dylib_plugin      = ./build_myrule/lib/libptlib_myrule.so   # path to the library
    phase_transitions = dylib                                   # builtin (default) | dylib
```

Both can also be given on the command line, which overrides the input file:
```
./bin/opt/LaMEM -ParamFile model.dat -dylib_plugin ./build_myrule/lib/libptlib_myrule.so -phase_transitions dylib
```

The library is loaded once at startup (also on a restart with `-mode restart`). LaMEM prints a parameter block for it, like for every other part of the setup, that states the library, its ABI version and whether the rule will actually be called:
```
Dylib plugin parameters:
   Library                                 : ./build_myrule/lib/libptlib_myrule.so
   Plugin ABI version                      : 2
   Phase transitions                       : dylib (lamem_phase_transition is called every step, after the built-in transitions)
--------------------------------------------------------------------------
```

With `phase_transitions = dylib` the rule is then called every time step, **after** any built-in `<PhaseTransitionStart>` blocks have been applied, and LaMEM reports what it did:
```
Dylib plugin  : 27512 marker(s) changed phase, 0 marker(s) changed temperature
```

A few points about how this interacts with the rest of the input file:

* The plugin does not require `Phasetrans = 1`. That flag only controls the built-in blocks; the plugin is controlled by `phase_transitions` alone. Built-in blocks and the plugin can be used together or separately.
* Because the plugin runs after the built-in transitions, it sees their results within the same time step, whereas the built-in transitions only see the plugin's changes on the following step. Keep this ordering in mind if a built-in block and the plugin act on the same phases.
* `phase_transitions = dylib` without a `dylib_plugin` is an error at startup. The reverse — a loaded plugin with `phase_transitions = builtin` — is allowed; the parameter block then shows `Phase transitions : builtin` and the rule is not called.

## Validation

The test `test/t40_PhaseTransitionPlugin` compares LaMEM runs that use the two example plugins against runs that use the equivalent built-in transitions; the residual histories must agree to the tolerances used for the other phase-transition tests. It is a good template for validating a new rule: implement it as a built-in transition where possible, then check that the plugin reproduces it before extending the rule beyond what the built-ins can do.

## Troubleshooting

**`phase_transitions = dylib requires a loaded plugin`** — add `dylib_plugin = <path>` to the input file or `-dylib_plugin <path>` on the command line.

**`dylib_plugin: could not open <path>`** or **`... not found next to <path>`** — the library or its bundle is incomplete. Check the path, and make sure the whole `build_<name>` directory was kept together.

**`dylib_plugin: plugin failed on at least one rank (rc=-1)`** — your rule threw an exception. Common causes are a method error for an unexpected value, indexing outside a lookup table, or attempting to print. Run the rule on a few `MarkerView` values in an ordinary Julia session to find it.

**`rc=-2` or `rc=-3`** — the plugin was built against a different version of `LaMEMPlugin.jl` than the one this LaMEM version expects. Rebuild it with the `LaMEMPlugin.jl` shipped alongside the LaMEM you are running.

**`non-finite T_out`** or **`out-of-range phase`** — the rule returned `NaN`/`Inf` for the temperature, or a phase ID that does not exist in the material table. LaMEM checks every marker and stops at the first bad one; the message includes its local index.

**The rule is never called** — check that `phase_transitions = dylib` is set; with the default `builtin`, LaMEM loads the plugin but the `Dylib plugin parameters` block shows `Phase transitions : builtin`.
