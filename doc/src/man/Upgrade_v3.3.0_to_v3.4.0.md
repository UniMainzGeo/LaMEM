# Upgrading from v3.3.0 to v3.4.0

LaMEM v3.4.0 is a **feature release with no breaking changes**. It adds user-defined phase
transitions written in Julia and loaded as a compiled plugin, and time-dependent (piecewise-constant)
velocities for internal velocity boxes. Both are opt-in. An unmodified v3.3.0 `.dat` file runs under
v3.4.0 and gives the same results. The release also contains two small bug fixes, and some changes to
the build system and CI that only matter if you compile LaMEM yourself.

---

## 0. TL;DR — the v3.4.0 upgrade checklist

**Your existing input files need no changes.** There are no renamed or removed keywords, no changes
to the default behaviour, and no changes to the restart format.

1. **New: Julia phase-transition plugins.** Write the rule in Julia, compile it once with JuliaC.jl
   into a shared library, and activate it with `dylib_plugin = <library>` and
   `phase_transitions = dylib`. Without these keywords nothing changes.
   **→ [§1](@ref "1. Julia phase-transition plugins")**
2. **New: time-dependent velocity boxes.** `vx`, `vy` and `vz` of a `<VelBoxStart>` block can now be
   given as a piecewise-constant function of time instead of a constant.
   **→ [§2](@ref "2. Time-dependent velocity boxes")**
3. **Fixed: two bugs.** One concerns the removal of merged markers (already in v3.3.1), the other the
   convergence warning of the constitutive-equation solver. **→ [§3](@ref "3. Bug fixes")**
4. **Building and testing.** Only relevant if you compile LaMEM yourself or run the test suite: the
   Makefile and the test scripts were extended to support Windows and macOS. Linux builds are
   unaffected. **→ [§4](@ref "4. Building, testing and CI")**
5. ⚠️ **Do not run a v3.4.0 input file that uses the new keywords with an older LaMEM.** The keywords
   do not exist there, so the older version has no code that reads them.
   **→ [§5](@ref "5. Pitfalls when upgrading to v3.4.0")**

The version number printed in the startup banner (and set in `Project.toml`) is now `3.4.0`.

---

## 1. Julia phase-transition plugins

The built-in phase transitions (`Constant`, `Clapeyron`, `Box`) cover the common cases. If you need a
rule they cannot express — one that depends on several marker quantities at once, a reaction with its
own kinetics, a lookup table — you can now write it in Julia. The rule is compiled ahead of time into
a shared library, which LaMEM loads at startup and calls once per time step on every MPI rank.
No Julia process is involved at run time.

The workflow has three steps. The complete guide, including the list of marker and cell fields that
a rule can read and write, is [Julia phase-transition plugins](JuliaPlugins.md).

**1. Write the rule** as a single Julia file that includes `scripts/dylib_plugins/LaMEMPlugin.jl`:

```julia
function my_rule(m::MarkerView)
    if m.phase == 2 || m.phase == 3
        return update(m; phase = m.T >= 1200.0 ? 3 : 2)
    end
    return m
end
```

**2. Build it** (needs Julia 1.12 or later and the `JuliaC` package):

```
julia --project=@juliac scripts/dylib_plugins/build_plugin.jl ptlib_myrule.jl
```

**3. Activate it** in the `.dat` file, or on the command line, which overrides the file:

```
    dylib_plugin      = ./build_myrule/lib/libptlib_myrule.so
    phase_transitions = dylib                     # builtin (default) | dylib
```

Points that are easy to miss:

* **Order within a time step.** The plugin runs *after* the built-in `<PhaseTransitionStart>` blocks.
  It therefore sees their result in the same step, whereas the built-in blocks only see the plugin's
  changes in the following step.
* **`Phasetrans = 1` is not needed.** That flag only controls the built-in blocks. The plugin is
  controlled by `phase_transitions` alone, and the two can be combined.
* **`phase_transitions = dylib` without `dylib_plugin` is an error** at startup. The reverse (a loaded
  plugin with `phase_transitions = builtin`) is allowed; the rule is then simply not called.
* **Restarts.** The library is loaded again on `-mode restart`; give the same `dylib_plugin`.
* **Limits.** One plugin per run. It works when LaMEM runs as its own process, including when it is
  started from `LaMEM.jl`'s `run_lamem`. It cannot be used from a LaMEM that is loaded as a library
  into a running Julia session.
* **Keep the whole bundle.** The library depends on the other files in its `build_<name>` directory.
  Copy the directory, never the library file alone.
* **Rebuild after updating.** The plugin interface is versioned (ABI v3). LaMEM refuses to load a
  plugin built against a different `LaMEMPlugin.jl` and says so at startup.

At startup LaMEM prints a `Dylib plugin parameters:` block with the library, its ABI version and
whether the rule is called, and once per step it reports how many markers changed phase, temperature
or other fields.

Platforms: Linux, macOS and Windows (tested in CI). The new keywords are listed in
`info/options/input_file.dat`.

---

## 2. Time-dependent velocity boxes

A `<VelBoxStart>` block prescribes a velocity inside a box. In v3.3.0 each of `vx`, `vy`, `vz` was a
constant. In v3.4.0 each component can instead be a piecewise-constant function of time, using the
same convention as the other `*_num_periods` keywords of the boundary conditions:

```
<VelBoxStart>
    cenX   =  4.0
    cenY   =  4.0
    cenZ   =  0.0
    widthX =  4.0
    widthY =  2.0
    widthZ =  2.0
    advect =  1

    vx_num_periods = 3            # number of periods
    vx_time_delims = 1.0 2.0      # end time of every period but the last (num_periods - 1 values)
    vx_values      = 1.0 0.0 -1.0 # velocity in each period
<VelBoxEnd>
```

Here `vx` is 1.0 until `t = 1.0`, then 0.0 until `t = 2.0`, then -1.0. The same pattern works for
`vy` (`vy_num_periods`, `vy_time_delims`, `vy_values`) and `vz`. The time delimiters use the time
unit and the values the velocity unit of the model (`units = geo`, `si` or `none`).

* **Per component.** One component can be time-dependent while another is constant, as in
  `vx_num_periods = 3` together with `vy = 0.4`.
* **Constant and periods are exclusive.** Giving both `vx` and `vx_num_periods` is an error. So are
  delimiters that are not in ascending order, and more than 20 periods.
* **`advect = 1`.** A box that moves with its velocity is displaced by the *integral* of the
  piecewise-constant velocity, not by `v*t`. For constant velocities the two are identical, so
  existing models are unchanged.
* **Backward compatible.** The block still has to specify at least one component, either a constant or
  periods. Existing `vx`/`vy`/`vz` input behaves exactly as in v3.3.0.

The test `test/t40_TimeDependentVBox` covers the fixed and the advected box, a box without periods,
and one with a time-dependent `vx` combined with a constant `vy`.

---

## 3. Bug fixes

* **Marker garbage collection (also in v3.3.1).** After marker resampling (cloning and merging of
  markers to keep the number of markers per cell within bounds), the markers that were merged away are
  removed from the marker storage. This requires their indices in ascending order, but the list was not
  sorted. It is sorted now. Runs with resampling and merging of markers are affected.
* **Convergence warning of the constitutive equations.** The `MPI_Reduce` that gathers the
  convergence statistics had its send and receive buffers swapped. The warning
  `Warning! Number of diverged points : ...` could therefore never trigger. It now does. **A run can
  newly print this warning in v3.4.0.** The warning is diagnostic only: the solution is not changed by
  the fix. If it appears, the local constitutive iterations did not converge at some integration
  points, and that was already the case in v3.3.0.

---

## 4. Building, testing and CI

None of the following changes the behaviour of a Linux build.

* **`src/Makefile`** now also builds when the C++ compiler is not an MPI wrapper (for example PETSc
  with Microsoft MPI on Windows), recognises GCC's `c++` front end, and handles MSYS2/MinGW (the
  linker flags and the PETSc libraries without an import library).
* **`src/LaMEM.h`, `src/Tensor.cpp`, `src/scaling.cpp`** define `__STDC_FORMAT_MACROS` and, on
  Windows, `_USE_MATH_DEFINES`, so that `PRIu64` and `M_PI` are available with every toolchain.
* **`src/dylib_plugins.{h,cpp}`** are new files. The Makefile picks them up without changes. Loading
  a plugin uses `PetscDLOpen`, so no new link dependency is added.
* **`test/compile_lamem.jl` and `test/setup_packages.jl`** can stage and compile against PETSc_jll on
  macOS and Windows, controlled by `LAMEM_CI_DESTDIR`. The default (`/workspace/destdir`) is unchanged.
* **CI** now also runs on macOS (arm64) and has a Windows job that builds LaMEM and runs one plugin
  test. The Linux job builds the plugin bundles that the new test needs. CI uses Julia 1.13.
* **New tests.** `t40_PhaseTransitionPlugin` compares plugin runs against the equivalent built-in
  transitions and checks that fields written by a plugin persist on the markers, also on two MPI
  ranks. It skips itself, with a message, if the plugin bundle has not been built. Run the
  plugin test locally by building the bundles first, as in the `Build dylib plugin bundle` step of
  `.github/workflows/CI.yml`.

---

## 5. Pitfalls when upgrading to v3.4.0

There are no silent changes to existing input files. The pitfalls concern the new features.

* **Running a v3.4.0 input file on an older version.** `phase_transitions`, `dylib_plugin` and the
  `vx_`/`vy_`/`vz_` period keywords do not exist before v3.4.0, so an older LaMEM has no code that
  reads them. A model that depends on a plugin or on a time-dependent box then does not run the
  physics you intended. Check the version in the startup banner.
* **Plugin built for another LaMEM.** After updating `LaMEMPlugin.jl`, rebuild every plugin. A
  mismatch is detected at load time (`reports ABI v2, expected v3`, `struct layout mismatch`) and
  stops the run with an error.
* **Plugin and built-in transitions on the same phases.** Because of the ordering in
  [§1](@ref "1. Julia phase-transition plugins"), a chain of transitions that would cascade within one
  step in the built-in code is delayed by one step when the second link is in the plugin.
* **Plugin in a Julia session.** A LaMEM that is loaded as a library into a running Julia session
  cannot load a plugin. Use `run_lamem`, which starts LaMEM as a subprocess.
