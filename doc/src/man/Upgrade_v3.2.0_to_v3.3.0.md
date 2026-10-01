# Upgrading from v3.2.0 to v3.3.0

LaMEM v3.3.0 is a **small feature and bug-fix release**. It adds mid-run phase injection for the
built-in geometric primitives (`n_inject` / `t_inject`), fixes the velocity on the bottom boundary
when the top is open, and brings the FastScape documentation in line with the code. There are no
build or test-harness changes. Most v3.2.0 `.dat` files run unchanged and give the same numbers.
The exceptions are open-boundary models with a background strain rate, and geometry blocks left
over in `msetup = files`/`polygons` inputs.

> **What this guide reflects.** The `v3.3.0` tag, which is also upstream `master` at commit
> `eb306fc4` (2026‑09‑29). No commits have landed after the tag. The baseline is `406b6444`
> (2026‑09‑22), the `v3.2.0` tag, which is the commit that *Upgrading from v3.1.0 to v3.2.0* was
> written against. The range is 14 commits and 4 merged PRs: #82 (phase injection), #88 (FastScape
> docs), #89 (the previous upgrade guide) and #90 (open-top bottom velocity), plus the version bump.
>
> Every claim was checked against both source trees. The behaviour in
> [§4](@ref "4. Pitfalls when upgrading to v3.3.0") was confirmed by running v3.2.0 and v3.3.0
> binaries side by side (PETSc 3.22.5, aarch64 Linux).

---

## 0. TL;DR — the v3.3.0 upgrade checklist

1. ⚠️ **Open-boundary models with a background strain rate give different results.**
   - **Open top:** with `open_top_bound = 1`, the bottom boundary now moves with the pure-shear
     velocity `vz = Ezz·(zbot − Rzz)`. In v3.0.0–v3.2.0 it was held at `vz = 0`.
   - **Open bottom:** with `open_bot_bound = 1` and a closed top, the bottom `vz` is now zero. In
     v3.2.0 it was the pure-shear value.
   - **Unchanged:** models with both boundaries open, or with `pres_bot` set.
   - Nothing is printed. **→ [§1](@ref "1. Open-boundary velocity fix in v3.3.0")**
2. **Restart databases do not carry over.** A v3.2.0 `restart/` directory cannot be read by
   v3.3.0, and a v3.3.0 one cannot be read by v3.2.0. The run aborts with a misleading PETSc
   ownership-range error. Finish or restart runs with the binary that wrote them.
   **→ [§4](@ref "4. Pitfalls when upgrading to v3.3.0")**
3. **`msetup = files` / `polygons`: geometry blocks are now parsed.** A leftover, complete block
   is still ignored, but now with a warning. An *incomplete or invalid* leftover block, for example
   one missing `bounds` or with an out-of-range `phase`, now stops the run. v3.2.0 silently ignored
   it. **→ [§2](@ref "2. Mid-run phase injection in v3.3.0")**
4. **New, opt-in: `n_inject` / `t_inject`** on every geometric primitive inject a body at given
   times during the run. **→ [§2](@ref "2. Mid-run phase injection in v3.3.0")**
5. **FastScape docs fixed.** `vel_boundary` is now documented correctly: `1` sets the velocity to
   zero and `0` keeps LaMEM's velocity. The code did not change. The warning in the v3.2.0 guide no
   longer applies to the docs. **→ [§3](@ref "3. FastScape documentation fixes in v3.3.0")**

### Quick self-check for a v3.2.0 input file

```bash
grep -nE "open_(top|bot)_bound\s*=\s*1" my_model.dat   # together with…
grep -nE "^\s*e(xx|yy)_(num_periods|strain_rates)" my_model.dat  # → both hit? results change (§1)
grep -nE "msetup\s*=\s*(files|polygons)" my_model.dat  # together with…
grep -nE "<(Sphere|Box|Layer|Cylinder|Ellipsoid|Hex|RidgeSeg)Start>" my_model.dat  # → leftover blocks now parsed (§2)
ls restart/ 2>/dev/null                                # → written by v3.2.0? cannot be restarted (§4)
```

No parameter was removed. LaMEM v3.2.0 could parse 517 parameters, and all 517 still parse in
v3.3.0. There are two new ones, `n_inject` and `t_inject` (`src/marker.cpp:871,875`).

---

## 1. Open-boundary velocity fix in v3.3.0

PR #90 changed `BCApplyVelDefault` (`src/bc.cpp:1432-1442`):

```cpp
// v3.2.0                                  // v3.3.0
if(top_open) { vez = 0.0; vbz = 0.0; }     if(top_open)     { vez = 0.0; }
                                           if(bc->bot_open) { vbz = 0.0; }
```

The mesh has always deformed with `Ezz` about `bg_ref_point`. In v3.2.0, an open top also held the
material at the bottom at `vz = 0`, so the mesh and the material came apart unless `bg_ref_point`
lay on the bottom boundary. Marker control then filled the gap. PR #90 reports that a basal salt
layer grew by 180 % in compression, and that the free surface subsided by about 20 km in a rift
model. The v3.3.0 behaviour matches the one for a closed top.

You are affected if you combine a background strain rate (`exx_*`/`eyy_*`, which set
`Ezz = −(Exx + Eyy)`) with exactly one open boundary: `open_top_bound = 1` alone, or
`open_bot_bound = 1` alone. The bottom must also not sit at the reference point (default
`bg_ref_point = 0 0 0`). With both boundaries open, `vbz` is zero in both versions. With `pres_bot`
set, the bottom `vz` is not constrained (`src/bc.cpp:1492`).

Measured in the first time step:

| Input | Version | bottom `vz` [cm/yr] | top `vz` [cm/yr] |
|---|---|---|---|
| `t14_1DStrengthEnvelope/1D_VP.dat`: open top, bottom at −54 km, `exx = −1e‑15` | v3.2.0 | 0 | 0.2019 |
| | v3.3.0 | −0.1704 | 0.0316 |
| `t23_Permeable/Permeable.dat` + `exx = −1e‑15`: open bottom, bottom at −1000 km | v3.2.0 | −3.156 | 0.158 |
| | v3.3.0 | 0 | 0.158 |

The same input with `open_top_bound = 0` gives byte-identical output in both versions. The t14
reference files and the four `norm(τII…)` targets in `test/runtests.jl` were updated. Results for
such models cannot be reproduced with v3.3.0. To get a fixed bottom, put `bg_ref_point` on
the bottom boundary, as PR #90 suggests; the mesh then deforms about the bottom as well.

---

## 2. Mid-run phase injection in v3.3.0

PR #82 adds two optional keys to every geometric primitive block. See the
[Built-in geometries](BuiltInGeometries.md) page and `info/options/input_file.dat`.

```
<SphereStart>
    phase    = 2
    center   = 20.0 50.0 80.0
    radius   = 15.0
    n_inject = 2          # number of injection times (default 0 = initial geometry)
    t_inject = 2.0 4.0    # inject at 2 and 4 Myr (GEO units)
<SphereEnd>
```

- **When it fires:** with `n_inject > 0`, the block is *not* applied at t = 0. At the start of the
  first step that reaches each time (`ADVMarkInjectGeom`, `src/marker.cpp:1378`, called from
  `LaMEMLib.cpp:637`), it overwrites the phase of the markers inside it. It also overwrites the
  temperature, if one is given.
- **Deformation history:** APS, ATS and the deviatoric stress of those markers are reset.
- **Limits:** up to 10 times (`_max_inj_times_`). They must be positive and strictly increasing,
  or the run stops with `t_inject values must be positive` or `… strictly increasing order`.
- **Restarts:** injections that have already fired are stored in the restart database.
- **Other setups:** it works for every `msetup`. For `files`/`polygons`, such as a
  GeophysicalModelGenerator setup, the marker file sets the initial geometry and primitive blocks
  are read only for their injections. Blocks without `n_inject` are ignored, with
  `Warning: N geometric primitive(s) without n_inject are ignored, …` (`src/marker.cpp:1369`).
  Because these blocks are now parsed, a leftover block that is incomplete or invalid fails, for
  example with `Define parameter "[-]bounds"`.
- **Test:** `t39_PhaseInjection` covers a full run and a restart. The next free test number is
  t40.

---

## 3. FastScape documentation fixes in v3.3.0

PR #88 changed only `doc/src/man/FastScape.md` and `info/options/input_file.dat`, not the code.
The FastScape pitfalls listed in *Upgrading from v3.1.0 to v3.2.0* are now documented where they
belong:
- `vel_boundary`: `1` sets the velocity to zero and `0` keeps it. The page also gives the digit
  order: bottom, right, top, left.
- `max_fs_dt` is in LaMEM time units.
- The FastScape output flags default to 1.
- `units = geo`/`si` and the `<FastScapeStart>` block are required.
- `surf_mode = 0` freezes the surface.
- The built-in surface keys are ignored with `surf_mode = 2`.

For source patchers: `AdvCtx` gained `numInjGeom` and `injGeom[_max_geom_]`
(`src/advect.h:144-145`), `GeomPrim` gained the injection fields (`src/marker.h`), and
`src/advect.h` now includes `marker.h`.

---

## 4. Pitfalls when upgrading to v3.3.0

### Silent: open boundaries with a background strain rate

The results change without any message. See [§1](@ref "1. Open-boundary velocity fix in v3.3.0").

### Silent: a v3.3.0 injection input run on v3.2.0

v3.2.0 ignores `n_inject`/`t_inject` without any message and places the bodies at t = 0. It also
accepts invalid values that v3.3.0 rejects.

### Silent: negative n_inject

`n_inject` is checked only against its maximum. `n_inject = -1` is accepted, and the primitive is
**never applied**, neither at t = 0 nor later. Use `0`, or leave the key out.

### Loud: restarting across versions

`LaMEMLib` is written to the restart file as one binary block (`src/LaMEMLib.cpp:252,334`), and
PR #82 made it larger. Reading a database written by the other version fails in
`FDSTAGReadRestart` with:

```
[0]PETSC ERROR: Ownership ranges sum to … but global dimension is …
```

There is no version check, so the message does not point to the real cause.

### Loud: leftover geometry blocks with msetup = files / polygons

A complete block now prints a warning and the output is unchanged. An incomplete or invalid block
now stops the run, for example with `Define parameter "[-]bounds"` or
`Entry 1 in parameter "[-]phase" is larger than allowed`. Delete the block, or fix it.

### Not a pitfall: the rest of your v3.2.0 inputs

Unmodified inputs without open boundaries, injection or leftover blocks gave byte-identical output
under both versions: t01 (falling block), t15 (RTI), t31 (polygons) and a closed-top t14 variant.
The full v3.3.0 suite (`make test` with `FASTSCAPE_LIB` set, so a `surf=scape` build) passes:
**107 passed, 0 failed** (25 min), including t37 in opt and deb mode and the new t39.
