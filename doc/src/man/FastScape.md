# Surface processes with FastScape

LaMEM can be coupled to the [FastScape](https://fastscape.org) (Fortran) library to compute landscape evolution on the internal free surface. Instead of the built-in erosion and sedimentation models, FastScape solves for fluvial incision, sediment transport and deposition, hillslope diffusion and — optionally — marine transport and compaction.

During each time step, LaMEM passes the surface velocities to FastScape, FastScape advances the topography (in one or more substeps), and the updated topography is returned to LaMEM, which then updates the free surface and converts markers between rock, air and sediment phases.

Both 2D and 3D LaMEM models are supported, on uniform and non-uniform grids.

!!! warning
    The FastScape surface grid is created once at startup and is never regenerated, so that drainage networks are not reset during the simulation. Setups in which the LaMEM grid stretches or contracts dynamically (i.e. models with a background strain rate) are therefore **currently not** supported.

!!! note
    FastScape runs serially on the first MPI rank, while LaMEM itself runs in parallel. This is usually cheap compared to the Stokes solve, but very large surface grids or large refinement factors may become a bottleneck.

## Compiling LaMEM with FastScape

FastScape support is not part of the default build; it has to be enabled explicitly.
This assumes you already have a working LaMEM installation with the PETSc version listed in the [Installation](Installation.md) section.

You need a compiled FastScape Fortran library (version 2.8.4 is what LaMEM is tested against). If you have Julia available, the easiest way to obtain a precompiled library for your platform is:
```julia
using Pkg
Pkg.add("Fastscapelib_jll")
using Fastscapelib_jll
println(Fastscapelib_jll.LIBPATH_list)
```
which prints the directory that contains `libfastscapelib_fortran`.

Set `FASTSCAPE_LIB` to the directory containing the library, in addition to the usual `PETSC_OPT`/`PETSC_DEB` variables:
```
export FASTSCAPE_LIB=/path/to/fastscape/lib
```

Then compile from the `/src` directory with the `surf=scape` option:
```
make mode=opt surf=scape all
```

This defines `-DWITH_FASTSCAPE` and links against `libfastscapelib_fortran`. You can check whether a given binary was built with FastScape support with:
```
./bin/opt/LaMEM -fastscape_info
```
which prints either `FASTSCAPE_ENABLED` or `FASTSCAPE_DISABLED`.

!!! note
    A binary compiled *without* FastScape support will silently fall back to the built-in surface processes if you set `surf_mode = 2`, rather than reporting an error. If you are unsure which binary you are using, check it with `-fastscape_info`.

## Activating FastScape in the input file

FastScape is activated by setting `surf_mode = 2` in the free-surface section. This requires the free surface itself to be active:
```
    surf_use       = 1     # free surface activation flag
    surf_level     = 0.0   # initial level
    surf_air_phase = 0     # phase ID of sticky air layer

    surf_mode      = 2     # 1-built-in erosion/sedimentation (default), 2-FastScape
```

Note that `surf_mode = 2` replaces the built-in surface processes: the `erosion_model` and `sediment_model` options are only used for `surf_mode = 1`.

All FastScape settings are then given in a separate block, usually placed near the end of the input file:
```
    <FastScapeStart>
        # --- Grid & coupling ---
        non_uniform_grid     = 0        # 0-uniform LaMEM grid, 1-non-uniform
        fs2D                 = 0        # 0-3D LaMEM model, 1-2D model
        fs_refine            = 2        # surface-grid refinement factor
        max_fs_dt            = 0.001    # maximum FastScape substep [Myr]
        sed_phases           = 3        # phase ID assigned to deposited sediment

        # --- Boundary conditions, initialization & output ---
        topo_boundary        = 1111     # topographic BC (4 digits)
        vel_boundary         = 0000     # surface velocity BC (4 digits)
        random_noise         = 1        # random perturbation of initial topography
        surf_out_nstep       = 1        # write FastScape output every n steps

        # --- Fluvial incision & hillslope transport ---
        kf                   = 1e-6     # bedrock river incision coefficient
        kfsed                = -1.0     # sediment incision coefficient (<0: use kf)
        m                    = 0.4      # drainage area exponent (stream power law)
        n                    = 1.0      # slope exponent (stream power law)
        kd                   = 1e-2     # bedrock hillslope diffusivity [m^2/yr]
        kdsed                = -1.0     # sediment diffusivity (<0: use kd)
        g                    = 0.0      # bedrock deposition coefficient
        gsed                 = 0.0      # sediment deposition coefficient (<0: use g)
        p                    = -2.0     # slope exponent for multi-direction flow

        set_marine           = 0        # marine transport & compaction [0-off, 1-on]
    <FastScapeEnd>
```

A complete, commented listing of **all** available parameters — including the marine transport options and the FastScape output flags — is given in the master input file [`info/options/input_file.dat`](https://github.com/UniMainzGeo/LaMEM/blob/master/info/options/input_file.dat).

A small working example is the test setup in [`test/t37_Collision_FastScape`](https://github.com/UniMainzGeo/LaMEM/tree/master/test/t37_Collision_FastScape), a continent-continent collision above a subducting slab, with the marker setup generated by `CreateMarkers_Collision.jl` using [GeophysicalModelGenerator.jl](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl).

### Key parameters

A few parameters deserve particular attention:

- `max_fs_dt` sets the maximum length of a FastScape substep. If the LaMEM time step is larger, it is subdivided into several FastScape substeps. Values that are too large under-resolve the landscape evolution values that are too small are simply expensive.
- `fs_refine = n` inserts `n-1` equally spaced points between adjacent surface nodes in each horizontal direction. Together with the LaMEM resolution, this controls whether the surface grid can resolve the expected valley spacing.
- `topo_boundary` and `vel_boundary` are four-digit codes, one digit per boundary. For the topography, `0` is a reflective (no-flux) boundary and `1` a fixed-height (open) boundary through which sediment leaves the domain. For the velocities, `0` sets the boundary velocity to zero and `1` keeps the velocity transferred from LaMEM.
- `sed_phases` is the LaMEM phase ID given to newly deposited sediment. The corresponding material properties have to be defined in the input file.

For a 2D LaMEM model (`fs2D = 1`), the 1D surface profile is extended in the second horizontal direction to build the 2D grid that FastScape requires. In that case `extendedRange` (and `extendedNodes`, for a uniform grid) must be given as well.

## Output

With `out_surf_fs = 1`, a coupled run writes an additional ParaView collection file next to the usual LaMEM output:

| File | Contents |
| :--- | :--- |
| `<out_file_name>.pvd` | main LaMEM (thermo-mechanical) output |
| `<out_file_name>_surf.pvd` | LaMEM free-surface output |
| `<out_file_name>_fs.pvd` | FastScape landscape output |

Which fields end up in the FastScape output is controlled by the `out_surf_*` flags, for example `out_surf_topofs`, `out_surf_erosion_rate`, `out_surf_drainage_area`, `out_surf_catchment` or `out_surf_lake_depth`. The write frequency is set with `surf_out_nstep`.

## Practical recommendations

- **Test the surface resolution.** The organisation of river networks can be sensitive to the FastScape grid resolution, which is a well-known limitation of most lanscape evolution models. Perform a resolution test before interpreting individual channels, drainage divides or small-scale landforms.
- **Test the coupling time step.** Vary `max_fs_dt` to check that the result is converged for your setup.
- **Use reproducible initial noise.** Drainage development depends on the initial topographic perturbation, so use the same `random_noise` setting when comparing runs with different physical parameters.
- **Choose boundary conditions deliberately.** Open (fixed-height) boundaries let sediment leave the domain, so mass is not conserved within the modelled surface; reflective boundaries retain sediment but are not appropriate everywhere.
- **Be careful with very regular drainage patterns.** Nearly uniform surface velocities and laterally homogeneous structures can produce unrealistically regular networks.

## Troubleshooting

**The build cannot find FastScape.**
Check that `FASTSCAPE_LIB` is set and visible in the current shell (`echo $FASTSCAPE_LIB`), and that it points to the directory that actually contains `libfastscapelib_fortran`. Also confirm that PETSc, FastScape and LaMEM were built with compatible compilers and MPI installations.

**FastScape does not seem to do anything.**
Verify that `surf_use = 1` and `surf_mode = 2`, and that the binary really has FastScape support (`LaMEM -fastscape_info`). A binary built without it falls back silently to the built-in surface processes.

**Sediment disappears through the boundaries.**
Review `topo_boundary`: open boundaries permit sediment discharge out of the domain.

**The coupled run is unexpectedly slow.**
Check the number of surface grid nodes (`fs_refine`, and `extendedNodes` in 2D) and the value of `max_fs_dt`. Both a highly refined surface grid and a small substep increase the serial FastScape workload.
