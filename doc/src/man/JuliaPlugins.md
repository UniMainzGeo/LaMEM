# User-defined phase transitions with Julia plugins

LaMEM's built-in phase transitions (`Constant`, `Clapeyron`, `Box`, see the `<PhaseTransitionStart>` blocks in `input_file.dat`) cover the common cases. When you need something they cannot express — a rule that depends on several marker quantities at once, a reaction with its own kinetics, a lookup table, a time-dependent threshold — you can write the rule in Julia instead and have LaMEM call it.

The rule is compiled once, ahead of time, with [JuliaC.jl](https://github.com/JuliaLang/JuliaC.jl) into a self-contained shared library (a *plugin*). LaMEM loads that library at startup and calls it once per time step, on every MPI rank, with all the markers that rank owns. There is no Julia process involved at run time and no interpreter overhead: the plugin is native code.

!!! note
    A plugin embeds its own Julia runtime. It works whenever LaMEM runs as its own process — including when it is launched from Julia through `LaMEM.jl`'s `run_lamem`, which starts LaMEM as a subprocess. It cannot be used from a LaMEM that is itself loaded as a library into an already-running Julia session. Only one plugin can be loaded per LaMEM run.

Plugins are supported on Linux, macOS and Windows, and are tested in CI on all three. You need Julia 1.12 or later and the `JuliaC` package (tested with 1.12 and 1.13).

## Writing a rule

A plugin is a single Julia file. It `include`s the helper module `scripts/dylib_plugins/LaMEMPlugin.jl` from the LaMEM repository, defines one function — the rule — and ends with a few lines of fixed boilerplate that export the rule to LaMEM.

The rule receives one marker at a time as a `MarkerView` and returns the marker, either unchanged or with some fields replaced using `update`:

```julia
module PTLibMyRule

include(joinpath(@__DIR__, "LaMEMPlugin.jl"))   # adjust the path to where you keep it
using .LaMEMPlugin

const PHASE_BELOW = Cint(2)
const PHASE_ABOVE = Cint(3)

# Called once per marker. Return `m` to leave the marker untouched, or
# `update(m; field = value, ...)` to change one or more writable fields.
function my_rule(m::MarkerView)
    if m.phase == PHASE_BELOW || m.phase == PHASE_ABOVE
        new_phase = m.T >= 1200.0 ? PHASE_ABOVE : PHASE_BELOW
        return update(m; phase = new_phase)
    end
    return m
end

# --- boilerplate, identical in every plugin (only the rule's name changes) ---
Base.@ccallable function lamem_phase_transition(markers::Ptr{LaMEMPluginMarkers}, cells::Ptr{LaMEMPluginCells},
        step::Ptr{LaMEMPluginStep}, scaling::Ptr{LaMEMPluginScaling})::Int32
    return lamem_pt_wrapper(my_rule, markers, cells, step, scaling)
end

end # module
```

`lamem_pt_wrapper` does all the bookkeeping: it converts LaMEM's internal non-dimensional data to dimensional values, builds a `MarkerView` for each marker, calls your rule, and writes back whatever the rule changed. Your rule only ever deals with physical quantities. `LaMEMPlugin.jl` also exports the two functions LaMEM uses to check, when it loads the library, that the plugin was built for the same interface (`lamem_plugin_abi_version` and `lamem_plugin_struct_sizes`); do not define them again in the rule file.

`update` takes any combination of the writable fields as keywords, for example `update(m; phase = 3, T = 900.0, aps = 0.0)`. For simple rules that only touch phase and temperature, the rule may instead return a tuple `(new_phase, new_T)`, as in plugins written for earlier LaMEM versions.

### What the rule sees

All values are **dimensional**, in the same units as your `.dat` file. The table lists them for `units = geo`; with `units = si` everything is SI, with `units = none` everything is non-dimensional.

**Marker fields** (`m.<field>`):

| Field | Meaning | Units (`geo`) | Writable |
|:--|:--|:--|:--|
| `phase` | phase ID (`Int32`) | | yes |
| `x`, `y`, `z` | position | km | no |
| `p` | pressure carried by the marker (the pressure history, see below; shifted like the pressure the rheology sees) | MPa | yes |
| `T` | temperature | °C | yes |
| `aps` | accumulated plastic strain | – | yes |
| `ats` | accumulated total strain | – | yes |
| `sxx`, `syy`, `szz`, `sxy`, `sxz`, `syz` | deviatoric (history) stress | MPa | yes |
| `ux`, `uy`, `uz` | total displacement | km | yes |
| `time`, `dt`, `step` | current time, time step size, step number | Myr | no |
| `T_internal` | input temperature in LaMEM's internal units, see below | – | no |
| `cell` | the cell the marker sits in, see below | | no |

**Cell fields** (`m.cell.<field>`): the values of the finite-difference cell the marker sits in, not interpolated to the marker position. They are read-only and hold the state at the start of the time step, i.e. the solution of the previous step (at the first step, solution-derived fields such as viscosity, stress and strain rate are not yet computed).

| Field | Meaning | Units (`geo`) |
|:--|:--|:--|
| `eta` | effective (total) viscosity | Pa·s |
| `eta_cr` | creep viscosity | Pa·s |
| `eta_st` | stabilization viscosity | Pa·s |
| `I2Gdt` | inverse elastic parameter 1/(2 G dt) | 1/(Pa·s) |
| `IKdt` | inverse bulk elastic parameter 1/(K dt) | 1/(Pa·s) |
| `Hr` | shear heating term | W/m³ |
| `Ha` | adiabatic heating term | W/m³ |
| `aps` | accumulated plastic strain | – |
| `ats` | accumulated total strain | – |
| `psr` | plastic strain-rate contribution (squared) | 1/s² |
| `theta` | volumetric strain rate | 1/s |
| `dxx`, `dyy`, `dzz` | total deviatoric strain rate | 1/s |
| `j2_strainrate` | second invariant of the deviatoric strain rate | 1/s |
| `sxx`, `syy`, `szz` | deviatoric stress | MPa |
| `hxx`, `hyy`, `hzz` | history (elastic) stress | MPa |
| `j2_stress` | second invariant of the deviatoric stress | MPa |
| `yield` | average yield stress | MPa |
| `rho` | density | kg/m³ |
| `rho_pf` | fluid density from a phase diagram | kg/m³ |
| `mf` | melt fraction | – |
| `alpha` | effective thermal expansion | 1/K |
| `cond` | thermal conductivity | W/m/K |
| `Tn` | temperature (history) | °C |
| `pn` | pressure (history) | MPa |
| `ux`, `uy`, `uz` | total displacement | km |
| `DIIdif`, `DIIdis`, `DIIprl`, `DIIfk`, `DIIpl` | relative diffusion, dislocation, Peierls, Frank-Kamenetzky and plastic strain rates | – |
| `phi` | PSD angle (adjoint only) | as stored |
| `free_surf` | free-surface flag (`Int32`) | – |
| `index` | the cell's local index on this MPI rank, 0-based (`Int32`) | – |

The phase composition of the cell is available as `phase_ratio(m.cell, ph)`, the volume fraction of phase `ph` (the phase ID as in the `.dat` file), and `num_phases(m.cell)`. Cell fields are read only when you access them, so a rule pays only for the fields it uses.

The `j2_*` invariants are cell-centred: the cell's own diagonal components plus the average of the four surrounding edge values for each off-diagonal component.

### What the rule may change

The **phase**, the **temperature**, the **pressure**, the **accumulated plastic and total strain** (`aps`, `ats`), the **deviatoric stress** (`sxx` ... `syz`) and the **displacement** (`ux`, `uy`, `uz`) of the marker. Only the position is read-only: moving a marker is the job of advection. Cell data is always read-only; change the markers, and LaMEM maps the changes back to the grid.

Writing the pressure needs some care, because the marker pressure is a *history* variable rather than the current solution. Each step, advection adds the change of the grid pressure to it; LaMEM then averages it onto the cells as the old pressure `pn` (`m.cell.pn`), which enters the next solve only through the volumetric elastic term of the continuity equation, `-(p - pn)/(K dt)`. Changing `p` therefore imposes a jump in the pressure history: for a compressible phase (bulk modulus `K` set) the next solve responds with volumetric expansion (lower `p`) or compaction (higher `p`); for an incompressible phase it has no effect. Density and the rheology use the current pressure, not this history value.

A changed field is converted back to internal units only if it actually changed, so any field your rule leaves alone is passed through bit-for-bit. Always start from `m` (`return m` or `update(m; ...)`) rather than recomputing an "unchanged" value, which could perturb it by a rounding error.

LaMEM rejects `NaN` or `Inf` in any written field and phase IDs outside the material table, and stops with an error naming the field and the marker.

### Constraints inside the rule

The plugin is compiled with `--trim=safe`, which keeps the library small but restricts what the rule may do:

* **No printing or other I/O** inside the rule. It is not supported in trimmed code and the plugin will fail to build.
* **Keep it type-stable.** Avoid dynamic dispatch, containers of abstract type, and global non-`const` state. Precompute constants with `const` at module level, as in the example. Access fields with literal names (`m.cell.eta`, not `getproperty(m.cell, name)` with a run-time `name`).
* **Exceptions abort the run.** The wrapper catches them and reports failure to LaMEM, which then stops with an error. Do not rely on exceptions for control flow.

### Reproducing a built-in transition exactly

LaMEM compares temperatures in its internal non-dimensional units. If you need a rule that reproduces a built-in transition bit-for-bit (for example to validate a plugin against a built-in run), compare `m.T_internal` against a threshold you non-dimensionalise with the same formula LaMEM uses, which `LaMEMPlugin` provides as `nondimensionalize_T`. `scripts/dylib_plugins/ptlib_constant.jl` does exactly this and documents the reasoning.

### The binary interface

You only need this section if you want to write a plugin in a language other than Julia, or change `LaMEMPlugin.jl`. The authoritative definition is `src/dylib_plugins.h` (ABI version 3). A plugin exports three C functions:

* `int32_t lamem_plugin_abi_version(void)` — must return 3.
* `int32_t lamem_plugin_struct_sizes(int64_t *sizes, int32_t n)` — writes the `sizeof` of the four structs below, in that order, and returns 4. LaMEM compares them with its own at load time.
* `int32_t lamem_phase_transition(const LaMEMPluginMarkers*, const LaMEMPluginCells*, const LaMEMPluginStep*, const LaMEMPluginScaling*)` — called once per step on every MPI rank, also when the rank owns no markers. It returns the number of markers it changed, or a negative value on failure.

The structs hold plain C pointers to arrays in LaMEM's internal units, together with the scaling factors to convert them (`LaMEMPluginScaling`). `LaMEMPluginMarkers` has one entry per marker of the rank: read-only `x`, `y`, `z`, the cell index `cell_index` (0-based), and an `*_in` / `*_out` array pair for every writable field (the pressure `p_in` / `p_out` without `pShift`). The `*_out` arrays are pre-filled with the input values, and LaMEM copies back only entries that differ from the input. `LaMEMPluginCells` holds one array per cell field for the rank's local cells, plus `phRat`, the phase ratios as one array of `ncells × numPhases` values (cell-major). `LaMEMPluginStep` holds the time, time step and step number.

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

The example rules in `scripts/dylib_plugins/` build the same way: `ptlib_constant.jl` and `ptlib_box.jl` reproduce the built-in `Constant` and `Box` transitions, and `ptlib_demo_fields.jl` shows a rule that reads cell data (temperature and phase ratio of the marker's cell) and writes fields other than phase and temperature (the accumulated plastic strain and the pressure).

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
   Plugin ABI version                      : 3
   Phase transitions                       : dylib (lamem_phase_transition is called every step, after the built-in transitions)
--------------------------------------------------------------------------
```

With `phase_transitions = dylib` the rule is then called every time step, **after** any built-in `<PhaseTransitionStart>` blocks have been applied, and LaMEM reports what it did:
```
Dylib plugin  : 27512 marker(s) changed phase, 0 marker(s) changed temperature, 0 marker(s) changed other fields
```

The last number counts markers on which the rule changed any of `p`, `aps`, `ats`, the stress or the displacement.

A few points about how this interacts with the rest of the input file:

* The plugin does not require `Phasetrans = 1`. That flag only controls the built-in blocks; the plugin is controlled by `phase_transitions` alone. Built-in blocks and the plugin can be used together or separately.
* Because the plugin runs after the built-in transitions, it sees their results within the same time step, whereas the built-in transitions only see the plugin's changes on the following step. Keep this ordering in mind if a built-in block and the plugin act on the same phases.
* `phase_transitions = dylib` without a `dylib_plugin` is an error at startup. The reverse — a loaded plugin with `phase_transitions = builtin` — is allowed; the parameter block then shows `Phase transitions : builtin` and the rule is not called.

## Validation

The test `test/t40_PhaseTransitionPlugin` compares LaMEM runs that use the `ptlib_constant` and `ptlib_box` plugins against runs that use the equivalent built-in transitions; the residual histories must agree to the tolerances used for the other phase-transition tests. It also runs `ptlib_demo_fields` on two MPI ranks and checks from the per-step summary that the plastic strain it writes persists on the markers, and checks the conversion of written values (including the pressure) back to internal units on hand-made input. It is a good template for validating a new rule: implement it as a built-in transition where possible, then check that the plugin reproduces it before extending the rule beyond what the built-ins can do.

## Troubleshooting

**`phase_transitions = dylib requires a loaded plugin`** — add `dylib_plugin = <path>` to the input file or `-dylib_plugin <path>` on the command line.

**`dylib_plugin: could not open <path>`** or **`... not found next to <path>`** — the library or its bundle is incomplete. Check the path, and make sure the whole `build_<name>` directory was kept together.

**`dylib_plugin: plugin failed on at least one rank (rc=-1)`** — your rule threw an exception. Common causes are a method error for an unexpected value, a phase ID outside the material table passed to `phase_ratio`, indexing outside a lookup table, or attempting to print. The rule can be tested in an ordinary Julia session by calling `lamem_pt_wrapper` on arrays you fill yourself.

**`reports ABI v2, expected v3`**, **`does not export lamem_plugin_abi_version`**, **`struct layout mismatch`**, **`rc=-2`** or **`rc=-3`** — the plugin was built against a different version of `LaMEMPlugin.jl` than the one this LaMEM version expects. Rebuild it with the `LaMEMPlugin.jl` shipped alongside the LaMEM you are running. Rules written for ABI v2 need their boilerplate replaced by the one shown above; a rule returning `(phase, T)` keeps working unchanged, but the `*_cell` fields of v2 are now `m.cell.j2_stress`, `m.cell.j2_strainrate`, `m.cell.eta` and `m.cell.aps`.

**`non-finite <field>_out`** or **`out-of-range phase`** — the rule wrote `NaN`/`Inf` into a field (for example `T_out` or `aps_out`), or a phase ID that does not exist in the material table. LaMEM checks every marker and stops at the first bad one; the message includes its local index.

**The rule is never called** — check that `phase_transitions = dylib` is set; with the default `builtin`, LaMEM loads the plugin but the `Dylib plugin parameters` block shows `Phase transitions : builtin`.
