# Shared helper for LaMEM phase-transition plugins (ABI v2, see
# src/phase_transition_plugin.h). A plugin `include`s this and calls
# lamem_pt_wrapper with a `rule(m::MarkerView) -> (phase, T_dimensional)`;
# all internal<->dimensional conversion happens here.
module LaMEMPlugin

export LaMEMPluginScaling, MarkerView, lamem_pt_wrapper,
       dimensionalize_T, nondimensionalize_T,
       dimensionalize_length, dimensionalize_time, dimensionalize_stress,
       dimensionalize_strain_rate, dimensionalize_viscosity,
       dimensionalize_velocity, dimensionalize_density,
       dimensionalize_pressure

const ABI_VERSION = Cint(2)

# mirrors src/phase_transition_plugin.h's LaMEMPluginScaling exactly
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
end

@inline dimensionalize_length(s::LaMEMPluginScaling, v)      = v * s.length
@inline dimensionalize_time(s::LaMEMPluginScaling, v)        = v * s.time
@inline dimensionalize_stress(s::LaMEMPluginScaling, v)      = v * s.stress
@inline dimensionalize_strain_rate(s::LaMEMPluginScaling, v) = v * s.strain_rate
@inline dimensionalize_viscosity(s::LaMEMPluginScaling, v)   = v * s.viscosity
@inline dimensionalize_velocity(s::LaMEMPluginScaling, v)    = v * s.velocity
@inline dimensionalize_density(s::LaMEMPluginScaling, v)     = v * s.density

@inline dimensionalize_T(s::LaMEMPluginScaling, T_internal) = T_internal * s.temperature - s.Tshift
@inline nondimensionalize_T(s::LaMEMPluginScaling, T_dim)   = (T_dim + s.Tshift) / s.temperature
@inline dimensionalize_pressure(s::LaMEMPluginScaling, p_internal_raw) = (p_internal_raw + s.pShift) * s.stress

# T_internal is kept alongside dimensional T so a plugin can reproduce a
# built-in comparison bit-for-bit in internal units (see ptlib_constant.jl)
struct MarkerView
    x::Float64
    y::Float64
    z::Float64
    T::Float64
    T_internal::Float64
    p::Float64
    time::Float64
    dt::Float64
    step::Int64
    sxx::Float64
    syy::Float64
    szz::Float64
    sxy::Float64
    sxz::Float64
    syz::Float64
    j2_stress_cell::Float64
    j2_strainrate_cell::Float64
    eta_cell::Float64
    aps_cell::Float64
    phase::Cint
end

# rule(m::MarkerView) -> (new_phase, new_T_dimensional), called once per marker
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

        # loud failure on a wrong struct layout, instead of silently reading garbage
        if s.abi_version != ABI_VERSION
            return Cint(-3)
        end
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
            Pout[i] = newph
            phase_changed = newph != Pin[i]

            # round-trip through dimensionalize/nondimensionalize only when
            # T actually changed: it is not bit-exact in general (Tshift is
            # not a power of 2), so an unconditional round-trip would
            # perturb P->T by ~1 ULP even when a rule never touches T
            T_changed = new_T_dim != m.T
            Tout[i] = T_changed ? nondimensionalize_T(s, Float64(new_T_dim)) : TT[i]

            (phase_changed || T_changed) && (changed += 1)
        end
        return Cint(changed)
    catch
        # uncaught exception in a @ccallable aborts the process; report -1 instead.
        # No I/O: printing is not --trim=safe inside a @ccallable function.
        return Cint(-1)
    end
end

lamem_phase_transition_abi_version() = ABI_VERSION

end # module
