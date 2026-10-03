module PTLibDemoFields

# Example rule that uses what the two re-implementations of built-in
# transitions (ptlib_constant.jl, ptlib_box.jl) do not: it reads cell data
# (m.cell) and writes a marker field other than phase or temperature.
#
# Rule: seed accumulated plastic strain (APS) in a hot plume head. A marker
# whose cell is hotter than T_SEED and consists of at least 50 % plume
# material (phase PLUME_PHASE) gets APS = APS_SEED, unless it already has at
# least that much. Seeded markers keep their APS when they move on, so after
# the first step only markers newly entering such a cell change. With a
# strain-softening rheology this would create a weak zone above the plume;
# in t40_PhaseTransitionPlugin (no plasticity) it only shows up in the APS
# output and in the "changed other fields" count LaMEM prints each step.
#
# The same markers also get their pressure lowered by P_DROP (a fraction),
# to exercise the pressure write path. The marker pressure is the pressure
# history (p_old): for a compressible phase this would act as a
# decompression in the next solve; the t40 model is incompressible, where
# it has no effect on the solution.
#
# A melt-weakening variant would look the same, e.g.
#     m.cell.mf > 0.1 ? update(m; aps = 0.0) : m
# (mf is only non-zero with phase diagrams or a melt parameterisation).

include(joinpath(@__DIR__, "LaMEMPlugin.jl"))
using .LaMEMPlugin

const T_SEED      = 1350.0 # cell temperature threshold (Celsius for units = geo)
const PLUME_PHASE = 4      # phase ID, as in the .dat file
const APS_SEED    = 1.0
const P_DROP      = 0.01   # relative pressure drop on seeded markers

function seed_rule(m::MarkerView)
    c = m.cell
    if m.aps < APS_SEED && c.Tn >= T_SEED && phase_ratio(c, PLUME_PHASE) >= 0.5
        return update(m; aps = APS_SEED, p = (1.0 - P_DROP) * m.p)
    end
    return m
end

Base.@ccallable function lamem_phase_transition(markers::Ptr{LaMEMPluginMarkers}, cells::Ptr{LaMEMPluginCells},
        step::Ptr{LaMEMPluginStep}, scaling::Ptr{LaMEMPluginScaling})::Int32
    return lamem_pt_wrapper(seed_rule, markers, cells, step, scaling)
end

end # module
