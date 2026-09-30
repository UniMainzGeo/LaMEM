module PTLibBox

# Reproduces the built-in Box transition used in box_builtin.dat: inside a
# box spanning the whole domain, phase 3 -> 2 (BothWays) and T is reset to a
# constant 900 C (Check_Box_Phase_Transition, PTBox_TempType=constant).

include("LaMEMPlugin.jl")
using .LaMEMPlugin

const XLO, XHI, YLO, YHI, ZLO, ZHI = -500.0, 500.0, -10.0, 10.0, -1000.0, 0.0
const PHASE_INSIDE, PHASE_OUTSIDE = Cint(2), Cint(3)
const CST_TEMP = 900.0

function box_rule(m::MarkerView)
    # Check_Box_Phase_Transition resets T for ANY marker geometrically
    # inside the box (not just phase 2/3); only phase 2/3 markers change phase.
    inside = XLO <= m.x <= XHI && YLO <= m.y <= YHI && ZLO <= m.z <= ZHI
    newT = inside ? CST_TEMP : m.T
    if m.phase == PHASE_INSIDE || m.phase == PHASE_OUTSIDE
        return (inside ? PHASE_INSIDE : PHASE_OUTSIDE, newT)
    end
    return (m.phase, newT)
end

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
        return lamem_pt_wrapper(box_rule, n, x, y, z, T, p, time,
            sxx, syy, szz, sxy, sxz, syz,
            j2_stress_cell, j2_strainrate_cell, eta_cell, aps_cell,
            phase_in, phase_out, T_out, scaling)
    catch
        return Cint(-1)
    end
end

Base.@ccallable function lamem_phase_transition_abi_version()::Cint
    return LaMEMPlugin.ABI_VERSION
end

end # module
