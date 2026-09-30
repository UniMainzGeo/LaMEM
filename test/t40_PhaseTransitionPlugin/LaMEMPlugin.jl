# Shared helper module for LaMEM phase-transition plugins (ABI v2). A user
# plugin `include`s this file and calls `lamem_pt_wrapper` from its own
# `Base.@ccallable lamem_phase_transition` function, passing a `rule`
# function that works entirely in DIMENSIONAL units - all internal <->
# dimensional unit conversion (using LaMEM's own characteristic scales,
# passed in by LaMEM every call) happens here, once, so a plugin author
# never has to touch scal->length/time/stress/... arithmetic themselves.
#
# See src/phase_transition_plugin.h (the authoritative ABI description) for
# the exact field layout/units of LaMEMPluginScaling and the exact
# conversion formulas this file implements.
module LaMEMPlugin

export LaMEMPluginScaling, MarkerView, lamem_pt_wrapper,
       dimensionalize_T, nondimensionalize_T,
       dimensionalize_length, dimensionalize_time, dimensionalize_stress,
       dimensionalize_strain_rate, dimensionalize_viscosity,
       dimensionalize_velocity, dimensionalize_density,
       dimensionalize_pressure

const ABI_VERSION = Cint(2)

# Mirrors src/phase_transition_plugin.h's `struct LaMEMPluginScaling`
# EXACTLY: same field order, same field types (isbits, C-compatible).
# int32_t -> Cint, double -> Cdouble, int64_t -> Clonglong (always 64-bit,
# unlike Clong which is 32-bit on Windows), char* -> Cstring/Ptr{Cchar}.
struct LaMEMPluginScaling
    abi_version::Cint
    utype::Cint
    length::Cdouble
    time::Cdouble
    stress::Cdouble
    temperature::Cdouble
    viscosity::Cdouble
    strain_rate::Cdouble
    velocity::Cdouble
    density::Cdouble
    Tshift::Cdouble
    pShift::Cdouble
    dt::Cdouble
    step::Clonglong
    lbl_length::Ptr{Cchar}
    lbl_time::Ptr{Cchar}
    lbl_stress::Ptr{Cchar}
    lbl_temperature::Ptr{Cchar}
    lbl_viscosity::Ptr{Cchar}
    lbl_strain_rate::Ptr{Cchar}
    lbl_velocity::Ptr{Cchar}
    lbl_density::Ptr{Cchar}
end

# --- unit conversion (internal <-> dimensional), see phase_transition_plugin.h ---
# All of these are the exact formulas LaMEM itself uses (src/scaling.h/.cpp,
# src/phase_transition.cpp's Set_Constant_Phase_Transition for the T one).

@inline dimensionalize_length(s::LaMEMPluginScaling, v)      = v * s.length
@inline dimensionalize_time(s::LaMEMPluginScaling, v)        = v * s.time
@inline dimensionalize_stress(s::LaMEMPluginScaling, v)      = v * s.stress
@inline dimensionalize_strain_rate(s::LaMEMPluginScaling, v) = v * s.strain_rate
@inline dimensionalize_viscosity(s::LaMEMPluginScaling, v)   = v * s.viscosity
@inline dimensionalize_velocity(s::LaMEMPluginScaling, v)    = v * s.velocity
@inline dimensionalize_density(s::LaMEMPluginScaling, v)     = v * s.density

# Temperature: dimensional = internal*temperature - Tshift (Celsius in geo
# mode, Kelvin in SI/none mode - Tshift is 0 there).
@inline dimensionalize_T(s::LaMEMPluginScaling, T_internal) = T_internal * s.temperature - s.Tshift
@inline nondimensionalize_T(s::LaMEMPluginScaling, T_dim)   = (T_dim + s.Tshift) / s.temperature

# Pressure the rheology/plasticity actually sees, dimensional (the raw `p`
# array element is P->p itself, WITHOUT pShift folded in - see
# phase_transition_plugin.h).
@inline dimensionalize_pressure(s::LaMEMPluginScaling, p_internal_raw) = (p_internal_raw + s.pShift) * s.stress

# --- per-marker dimensional view, built once per marker inside the wrapper loop ---
# A plain isbits struct (not a closure/generic wrapper) so it stays
# trim-safe: every field is a concrete Float64/Int32, no Any-typed fields.
struct MarkerView
    x::Float64
    y::Float64
    z::Float64
    T::Float64              # dimensional (Celsius in geo mode)
    T_internal::Float64     # LaMEM's raw internal (non-dimensional) T, i.e.
                             # P->T itself, unconverted - kept alongside the
                             # dimensional T so a plugin that wants to
                             # reproduce a built-in comparison bit-for-bit in
                             # INTERNAL units (see ptlib_constant.jl) can,
                             # without having to invert `T` back itself
    p::Float64              # dimensional, rheology-seen pressure (pShift folded in)
    time::Float64           # dimensional simulation time
    dt::Float64             # dimensional time step
    step::Int64
    sxx::Float64
    syy::Float64
    szz::Float64
    sxy::Float64
    sxz::Float64
    syz::Float64
    j2_stress_cell::Float64      # dimensional
    j2_strainrate_cell::Float64  # dimensional
    eta_cell::Float64            # dimensional
    aps_cell::Float64            # dimensionless
    phase::Cint
end

# The wrapper: LaMEM calls this (via the plugin's own
# Base.@ccallable lamem_phase_transition, which simply forwards its
# arguments here) once per time step, once per MPI rank, with n local
# markers. `rule(m::MarkerView) -> (new_phase::Integer, new_T_dimensional::Real)`
# is called once per marker; return `(m.phase, m.T)` unchanged for a marker
# the rule does not want to touch.
function lamem_pt_wrapper(rule::F,
        n::Csize_t,
        x::Ptr{Cdouble}, y::Ptr{Cdouble}, z::Ptr{Cdouble},
        T::Ptr{Cdouble}, p::Ptr{Cdouble}, time::Cdouble,
        sxx::Ptr{Cdouble}, syy::Ptr{Cdouble}, szz::Ptr{Cdouble},
        sxy::Ptr{Cdouble}, sxz::Ptr{Cdouble}, syz::Ptr{Cdouble},
        j2_stress_cell::Ptr{Cdouble}, j2_strainrate_cell::Ptr{Cdouble},
        eta_cell::Ptr{Cdouble}, aps_cell::Ptr{Cdouble},
        phase_in::Ptr{Cint}, phase_out::Ptr{Cint},
        T_out::Ptr{Cdouble},
        scaling::Ptr{LaMEMPluginScaling})::Cint where F
    try
        nn = Int(n)
        s = unsafe_load(scaling)

        # Sanity check on the scaling struct itself: utype/length/time/stress
        # must always be positive for any real LaMEM run (length/time/stress
        # are characteristic SCALES, never zero or negative; utype is one of
        # 0/1/2 - see phase_transition_plugin.h - so utype<0 is impossible
        # for a genuinely valid struct too, though utype itself can
        # legitimately be 0 for _NONE_, so only length/time/stress are
        # checked for strict positivity here). This exists so that a wrong
        # struct layout (e.g. a future libblastrampoline-style ABI mismatch,
        # a stale plugin built against a different field order) fails LOUDLY
        # with a distinctive return code (-2) instead of silently reading
        # garbage and producing wrong-but-plausible-looking results.
        if s.utype < 0 || !(s.length > 0.0) || !(s.time > 0.0) || !(s.stress > 0.0)
            return Cint(-2)
        end

        X   = unsafe_wrap(Array, x, nn); Y = unsafe_wrap(Array, y, nn); Z = unsafe_wrap(Array, z, nn)
        TT  = unsafe_wrap(Array, T, nn); PP = unsafe_wrap(Array, p, nn)
        SXX = unsafe_wrap(Array, sxx, nn); SYY = unsafe_wrap(Array, syy, nn); SZZ = unsafe_wrap(Array, szz, nn)
        SXY = unsafe_wrap(Array, sxy, nn); SXZ = unsafe_wrap(Array, sxz, nn); SYZ = unsafe_wrap(Array, syz, nn)
        J2S = unsafe_wrap(Array, j2_stress_cell, nn); J2E = unsafe_wrap(Array, j2_strainrate_cell, nn)
        ETA = unsafe_wrap(Array, eta_cell, nn); APS = unsafe_wrap(Array, aps_cell, nn)
        Pin  = unsafe_wrap(Array, phase_in, nn)
        Pout = unsafe_wrap(Array, phase_out, nn)
        Tout = unsafe_wrap(Array, T_out, nn)

        time_dim = dimensionalize_time(s, time)
        dt_dim   = dimensionalize_time(s, s.dt)

        changed = 0
        @inbounds for i in 1:nn
            m = MarkerView(
                dimensionalize_length(s, X[i]),
                dimensionalize_length(s, Y[i]),
                dimensionalize_length(s, Z[i]),
                dimensionalize_T(s, TT[i]),
                TT[i],
                dimensionalize_pressure(s, PP[i]),
                time_dim,
                dt_dim,
                Int64(s.step),
                dimensionalize_stress(s, SXX[i]),
                dimensionalize_stress(s, SYY[i]),
                dimensionalize_stress(s, SZZ[i]),
                dimensionalize_stress(s, SXY[i]),
                dimensionalize_stress(s, SXZ[i]),
                dimensionalize_stress(s, SYZ[i]),
                dimensionalize_stress(s, J2S[i]),
                dimensionalize_strain_rate(s, J2E[i]),
                dimensionalize_viscosity(s, ETA[i]),
                APS[i],
                Pin[i],
            )

            new_phase, new_T_dim = rule(m)

            newph = Cint(new_phase)
            if newph != Pin[i]
                changed += 1
            end
            Pout[i] = newph

            # IMPORTANT: only round-trip T through dimensionalize/
            # nondimensionalize when the rule actually changed it (compared
            # in DIMENSIONAL units, against m.T). (internal -> dimensional
            # -> internal) is NOT bit-exact in general for an arbitrary
            # internal value (Tshift is not a power of 2, so the
            # intermediate *temperature - Tshift / +Tshift /temperature
            # arithmetic rounds differently for a large fraction of inputs -
            # empirically about 1 in 4 uniformly-random values fails to
            # round-trip exactly). Writing the round-tripped value back
            # UNCONDITIONALLY would silently perturb P->T by ~1 ULP on a
            # large fraction of markers every single time step, even for a
            # rule that never touches T at all - passing the ORIGINAL
            # internal value straight through instead keeps a "do nothing"
            # rule bit-for-bit inert, which is essential for the t40 test
            # (which asserts the built-in and plugin runs agree exactly).
            if new_T_dim == m.T
                Tout[i] = TT[i]
            else
                Tout[i] = nondimensionalize_T(s, Float64(new_T_dim))
            end
        end
        return Cint(changed)
    catch
        # An uncaught exception inside a @ccallable function aborts the
        # whole process outside PETSc's error handling, so every exception
        # is caught here and turned into a negative return value instead
        # (LaMEM turns a negative return into a normal, collective SETERRQ -
        # see phase_transition_plugin.h). No I/O (e.g. printing the
        # exception) happens here: even a fixed-string println pulls in
        # dynamic dispatch the --trim=safe verifier cannot resolve inside a
        # @ccallable function, so the build would fail.
        return Cint(-1)
    end
end

# Optional ABI version symbol - see phase_transition_plugin.h's "ABI
# VERSIONING". A plugin that `include`s this module and re-exports this
# function as its own @ccallable gets the version check for free.
lamem_phase_transition_abi_version() = ABI_VERSION

end # module
