# Shared helper for LaMEM phase-transition plugins (ABI v3, see
# src/dylib_plugins.h). A plugin `include`s this and calls
# `lamem_pt_wrapper(rule, markers, cells, step, scaling)` from its
# `lamem_phase_transition` @ccallable, with `rule(m::MarkerView)` returning
# the (possibly modified) marker; all internal<->dimensional conversion
# happens here. This module also exports the two other symbols LaMEM
# requires, lamem_plugin_abi_version and lamem_plugin_struct_sizes.
module LaMEMPlugin

export LaMEMPluginMarkers, LaMEMPluginCells, LaMEMPluginStep, LaMEMPluginScaling,
       MarkerView, CellView, lamem_pt_wrapper, update, phase_ratio, num_phases,
       dimensionalize_T, nondimensionalize_T,
       dimensionalize_length, dimensionalize_time, dimensionalize_stress,
       dimensionalize_strain_rate, dimensionalize_viscosity,
       dimensionalize_velocity, dimensionalize_density,
       dimensionalize_pressure

const ABI_VERSION = Int32(3)

#-----------------------------------------------------------------------------
# ABI structs: mirror src/dylib_plugins.h field-by-field (same order, same
# types). LaMEM compares their sizeof()s with its own at load time, see
# lamem_plugin_struct_sizes below.
#-----------------------------------------------------------------------------
struct LaMEMPluginMarkers
    n::Csize_t
    cell_index::Ptr{Int32}

    # read-only
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    z::Ptr{Cdouble}
    p::Ptr{Cdouble}

    # writable fields: input values
    phase_in::Ptr{Int32}
    T_in::Ptr{Cdouble}
    aps_in::Ptr{Cdouble}
    ats_in::Ptr{Cdouble}
    sxx_in::Ptr{Cdouble}
    syy_in::Ptr{Cdouble}
    szz_in::Ptr{Cdouble}
    sxy_in::Ptr{Cdouble}
    sxz_in::Ptr{Cdouble}
    syz_in::Ptr{Cdouble}
    ux_in::Ptr{Cdouble}
    uy_in::Ptr{Cdouble}
    uz_in::Ptr{Cdouble}

    # writable fields: output values (pre-filled with the input values)
    phase_out::Ptr{Int32}
    T_out::Ptr{Cdouble}
    aps_out::Ptr{Cdouble}
    ats_out::Ptr{Cdouble}
    sxx_out::Ptr{Cdouble}
    syy_out::Ptr{Cdouble}
    szz_out::Ptr{Cdouble}
    sxy_out::Ptr{Cdouble}
    sxz_out::Ptr{Cdouble}
    syz_out::Ptr{Cdouble}
    ux_out::Ptr{Cdouble}
    uy_out::Ptr{Cdouble}
    uz_out::Ptr{Cdouble}
end

struct LaMEMPluginCells
    ncells::Csize_t
    numPhases::Int32
    reserved::Int32

    # SolVarCell::svDev
    eta::Ptr{Cdouble}
    eta_st::Ptr{Cdouble}
    I2Gdt::Ptr{Cdouble}
    Hr::Ptr{Cdouble}
    aps::Ptr{Cdouble}
    psr::Ptr{Cdouble}

    # SolVarCell::svBulk
    theta::Ptr{Cdouble}
    rho::Ptr{Cdouble}
    IKdt::Ptr{Cdouble}
    alpha::Ptr{Cdouble}
    Tn::Ptr{Cdouble}
    pn::Ptr{Cdouble}
    rho_pf::Ptr{Cdouble}
    mf::Ptr{Cdouble}
    phi::Ptr{Cdouble}
    Ha::Ptr{Cdouble}
    cond::Ptr{Cdouble}

    # remaining SolVarCell fields
    sxx::Ptr{Cdouble}
    syy::Ptr{Cdouble}
    szz::Ptr{Cdouble}
    hxx::Ptr{Cdouble}
    hyy::Ptr{Cdouble}
    hzz::Ptr{Cdouble}
    dxx::Ptr{Cdouble}
    dyy::Ptr{Cdouble}
    dzz::Ptr{Cdouble}
    free_surf::Ptr{Int32}
    ux::Ptr{Cdouble}
    uy::Ptr{Cdouble}
    uz::Ptr{Cdouble}
    ats::Ptr{Cdouble}
    eta_cr::Ptr{Cdouble}
    DIIdif::Ptr{Cdouble}
    DIIdis::Ptr{Cdouble}
    DIIprl::Ptr{Cdouble}
    DIIfk::Ptr{Cdouble}
    DIIpl::Ptr{Cdouble}
    yield::Ptr{Cdouble}

    # derived: cell-centred J2 invariants
    j2_stress::Ptr{Cdouble}
    j2_strainrate::Ptr{Cdouble}

    # ncells*numPhases, cell-major
    phRat::Ptr{Cdouble}
end

struct LaMEMPluginStep
    time::Cdouble
    dt::Cdouble
    step::Int64
end

struct LaMEMPluginScaling
    abi_version::Int32
    utype::Int32
    length::Cdouble
    time::Cdouble
    stress::Cdouble
    temperature::Cdouble
    viscosity::Cdouble
    strain_rate::Cdouble
    velocity::Cdouble
    density::Cdouble
    conductivity::Cdouble
    expansivity::Cdouble
    dissipation_rate::Cdouble
    Tshift::Cdouble
    pShift::Cdouble
end

#-----------------------------------------------------------------------------
# unit conversion (internal -> dimensional, in the .dat file's units)
#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
# CellView: lazy, dimensional, read-only view of one cell
#-----------------------------------------------------------------------------
# Holds only pointers and the 1-based cell index, so building one per marker
# is free; a field is read (and dimensionalised) when it is accessed.
struct CellView
    _cells::Ptr{LaMEMPluginCells}
    _scaling::Ptr{LaMEMPluginScaling}
    _i::Int
end

const CELL_FIELDS = (:eta, :eta_st, :I2Gdt, :Hr, :aps, :psr,
    :theta, :rho, :IKdt, :alpha, :Tn, :pn, :rho_pf, :mf, :phi, :Ha, :cond,
    :sxx, :syy, :szz, :hxx, :hyy, :hzz, :dxx, :dyy, :dzz,
    :free_surf, :ux, :uy, :uz, :ats, :eta_cr,
    :DIIdif, :DIIdis, :DIIprl, :DIIfk, :DIIpl, :yield,
    :j2_stress, :j2_strainrate, :index)

Base.propertynames(::CellView, ::Bool=false) = CELL_FIELDS

# internal -> dimensional for a cell field (unit table: doc/src/man/JuliaPlugins.md)
@inline function _cell_dim(s::LaMEMPluginScaling, name::Symbol, v::Float64)
    if name === :eta || name === :eta_st || name === :eta_cr
        return v * s.viscosity
    elseif name === :I2Gdt || name === :IKdt
        return v / s.viscosity
    elseif name === :Hr || name === :Ha
        return v * s.dissipation_rate
    elseif name === :psr
        return v * s.strain_rate * s.strain_rate
    elseif name === :theta || name === :dxx || name === :dyy || name === :dzz || name === :j2_strainrate
        return v * s.strain_rate
    elseif name === :rho || name === :rho_pf
        return v * s.density
    elseif name === :alpha
        return v * s.expansivity
    elseif name === :Tn
        return dimensionalize_T(s, v)
    elseif name === :pn
        return dimensionalize_pressure(s, v)
    elseif name === :cond
        return v * s.conductivity
    elseif name === :sxx || name === :syy || name === :szz ||
           name === :hxx || name === :hyy || name === :hzz ||
           name === :yield || name === :j2_stress
        return v * s.stress
    elseif name === :ux || name === :uy || name === :uz
        return v * s.length
    else # aps, ats, mf, phi, DIIdif, DIIdis, DIIprl, DIIfk, DIIpl: no units
        return v
    end
end

Base.@constprop :aggressive @inline function Base.getproperty(c::CellView, name::Symbol)
    i = getfield(c, :_i)
    name === :index && return Int32(i - 1) # 0-based, as LaMEM's cell_index
    cells = unsafe_load(getfield(c, :_cells))
    name === :free_surf && return unsafe_load(cells.free_surf, i)
    v = unsafe_load(getfield(cells, name)::Ptr{Cdouble}, i)
    return _cell_dim(unsafe_load(getfield(c, :_scaling)), name, v)
end

# volume fraction of phase `ph` (0-based phase ID, as in the .dat file)
@inline function phase_ratio(c::CellView, ph::Integer)
    cells = unsafe_load(getfield(c, :_cells))
    np = Int(cells.numPhases)
    (0 <= ph < np) || throw(BoundsError())
    return unsafe_load(cells.phRat, (getfield(c, :_i) - 1) * np + Int(ph) + 1)
end

@inline num_phases(c::CellView) = Int(unsafe_load(getfield(c, :_cells)).numPhases)

#-----------------------------------------------------------------------------
# MarkerView: one marker, all fields dimensional
#-----------------------------------------------------------------------------
# T_internal is the input temperature in LaMEM's internal units, so a plugin
# can reproduce a built-in comparison bit-for-bit (see ptlib_constant.jl).
# Writable: phase, T, aps, ats, sxx..syz, ux..uz (return them changed via
# `update`); everything else is ignored on return.
struct MarkerView
    x::Float64
    y::Float64
    z::Float64
    p::Float64
    T::Float64
    aps::Float64
    ats::Float64
    sxx::Float64
    syy::Float64
    szz::Float64
    sxy::Float64
    sxz::Float64
    syz::Float64
    ux::Float64
    uy::Float64
    uz::Float64
    phase::Int32
    T_internal::Float64
    time::Float64
    dt::Float64
    step::Int64
    cell::CellView
end

# copy of `m` with the given writable fields replaced, e.g.
#   update(m; phase = 3, T = 900.0)
@inline function update(m::MarkerView;
        phase::Integer = m.phase, T::Real = m.T, aps::Real = m.aps, ats::Real = m.ats,
        sxx::Real = m.sxx, syy::Real = m.syy, szz::Real = m.szz,
        sxy::Real = m.sxy, sxz::Real = m.sxz, syz::Real = m.syz,
        ux::Real = m.ux, uy::Real = m.uy, uz::Real = m.uz)
    return MarkerView(m.x, m.y, m.z, m.p,
        Float64(T), Float64(aps), Float64(ats),
        Float64(sxx), Float64(syy), Float64(szz), Float64(sxy), Float64(sxz), Float64(syz),
        Float64(ux), Float64(uy), Float64(uz),
        Int32(phase), m.T_internal, m.time, m.dt, m.step, m.cell)
end

# a rule may also return the v2-style tuple (new_phase, new_T)
@inline _as_marker(r::MarkerView, ::MarkerView) = r
@inline _as_marker(r::Tuple{Integer, Real}, m::MarkerView) = update(m; phase = r[1], T = r[2])

# store `new` (dimensional) converted back to internal units, but only if the
# rule changed it: the round-trip is not bit-exact in general (e.g. Tshift is
# not a power of 2), and an unchanged entry must stay exactly as LaMEM
# pre-filled it. NaN/Inf are stored and rejected by LaMEM.
@inline function _store!(out::Ptr{Cdouble}, i::Int, new::Float64, old::Float64, scale::Float64)
    new == old && return false
    unsafe_store!(out, new / scale, i)
    return true
end

@inline function _store_T!(out::Ptr{Cdouble}, i::Int, new::Float64, old::Float64, s::LaMEMPluginScaling)
    new == old && return false
    unsafe_store!(out, nondimensionalize_T(s, new), i)
    return true
end

#-----------------------------------------------------------------------------
# entry point
#-----------------------------------------------------------------------------
# rule(m::MarkerView) -> MarkerView (or (phase, T)), called once per marker.
# Returns the number of markers changed, or -1 (exception), -2 (implausible
# scaling struct), -3 (ABI mismatch).
function lamem_pt_wrapper(rule::F,
        markers::Ptr{LaMEMPluginMarkers}, cells::Ptr{LaMEMPluginCells},
        step::Ptr{LaMEMPluginStep}, scaling::Ptr{LaMEMPluginScaling})::Int32 where F
    try
        s = unsafe_load(scaling)

        # loud failure on a wrong struct layout, instead of silently reading garbage
        if s.abi_version != ABI_VERSION
            return Int32(-3)
        end
        if s.utype < 0 || !(s.length > 0.0) || !(s.time > 0.0) || !(s.stress > 0.0) || !(s.temperature > 0.0)
            return Int32(-2)
        end

        M  = unsafe_load(markers)
        st = unsafe_load(step)
        n  = Int(M.n)

        time_dim = dimensionalize_time(s, st.time)
        dt_dim   = dimensionalize_time(s, st.dt)

        changed = 0
        @inbounds for i in 1:n
            T_int = unsafe_load(M.T_in, i)
            m = MarkerView(
                dimensionalize_length(s, unsafe_load(M.x, i)),
                dimensionalize_length(s, unsafe_load(M.y, i)),
                dimensionalize_length(s, unsafe_load(M.z, i)),
                dimensionalize_pressure(s, unsafe_load(M.p, i)),
                dimensionalize_T(s, T_int),
                unsafe_load(M.aps_in, i),
                unsafe_load(M.ats_in, i),
                dimensionalize_stress(s, unsafe_load(M.sxx_in, i)),
                dimensionalize_stress(s, unsafe_load(M.syy_in, i)),
                dimensionalize_stress(s, unsafe_load(M.szz_in, i)),
                dimensionalize_stress(s, unsafe_load(M.sxy_in, i)),
                dimensionalize_stress(s, unsafe_load(M.sxz_in, i)),
                dimensionalize_stress(s, unsafe_load(M.syz_in, i)),
                dimensionalize_length(s, unsafe_load(M.ux_in, i)),
                dimensionalize_length(s, unsafe_load(M.uy_in, i)),
                dimensionalize_length(s, unsafe_load(M.uz_in, i)),
                unsafe_load(M.phase_in, i),
                T_int,
                time_dim,
                dt_dim,
                st.step,
                CellView(cells, scaling, Int(unsafe_load(M.cell_index, i)) + 1),
            )

            r = _as_marker(rule(m), m)

            ch = false
            if r.phase != m.phase
                unsafe_store!(M.phase_out, r.phase, i)
                ch = true
            end
            ch |= _store_T!(M.T_out, i, r.T, m.T, s)
            ch |= _store!(M.aps_out, i, r.aps, m.aps, 1.0)
            ch |= _store!(M.ats_out, i, r.ats, m.ats, 1.0)
            ch |= _store!(M.sxx_out, i, r.sxx, m.sxx, s.stress)
            ch |= _store!(M.syy_out, i, r.syy, m.syy, s.stress)
            ch |= _store!(M.szz_out, i, r.szz, m.szz, s.stress)
            ch |= _store!(M.sxy_out, i, r.sxy, m.sxy, s.stress)
            ch |= _store!(M.sxz_out, i, r.sxz, m.sxz, s.stress)
            ch |= _store!(M.syz_out, i, r.syz, m.syz, s.stress)
            ch |= _store!(M.ux_out,  i, r.ux,  m.ux,  s.length)
            ch |= _store!(M.uy_out,  i, r.uy,  m.uy,  s.length)
            ch |= _store!(M.uz_out,  i, r.uz,  m.uz,  s.length)

            changed += ch
        end
        return Int32(changed)
    catch
        # uncaught exception in a @ccallable aborts the process; report -1 instead.
        # No I/O: printing is not --trim=safe inside a @ccallable function.
        return Int32(-1)
    end
end

#-----------------------------------------------------------------------------
# symbols LaMEM checks at load time (exported from every plugin that
# includes this file; a rule file must not define them again)
#-----------------------------------------------------------------------------
Base.@ccallable function lamem_plugin_abi_version()::Int32
    return ABI_VERSION
end

Base.@ccallable function lamem_plugin_struct_sizes(sizes::Ptr{Int64}, n::Int32)::Int32
    sz = (sizeof(LaMEMPluginMarkers), sizeof(LaMEMPluginCells),
          sizeof(LaMEMPluginStep), sizeof(LaMEMPluginScaling))
    for i in 1:min(Int(n), length(sz))
        unsafe_store!(sizes, Int64(sz[i]), i)
    end
    return Int32(length(sz))
end

end # module
