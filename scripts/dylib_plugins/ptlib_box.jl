module PTLibBox

# Reproduces the built-in Box transition used in box_builtin.dat: inside a
# box spanning the whole domain, phase 3 -> 2 (BothWays) and T is reset to a
# constant 900 C (Check_Box_Phase_Transition, PTBox_TempType=constant).

include(joinpath(@__DIR__, "LaMEMPlugin.jl"))
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
        return update(m; phase = inside ? PHASE_INSIDE : PHASE_OUTSIDE, T = newT)
    end
    return update(m; T = newT)
end

Base.@ccallable function lamem_phase_transition(markers::Ptr{LaMEMPluginMarkers}, cells::Ptr{LaMEMPluginCells},
        step::Ptr{LaMEMPluginStep}, scaling::Ptr{LaMEMPluginScaling})::Int32
    return lamem_pt_wrapper(box_rule, markers, cells, step, scaling)
end

end # module
