module PTLibConstant

# Julia re-implementation of LaMEM's built-in "Constant" phase transition,
# as configured by PhaseTransition ID 0 in
# test/t16_PhaseTransitions/Plume_PhaseTransitions.dat:
#
#   Type                  = Constant
#   Parameter_transition  = T
#   ConstantValue         = 1200          # dimensional, scal units (Celsius in geo mode)
#   PhaseAbove            = 3
#   PhaseBelow            = 2
#   PhaseDirection        = BothWays
#
# LaMEM's Check_Constant_Phase_Transition (src/phase_transition.cpp) rule for
# Parameter_transition == T is, for a marker whose phase is PhaseBelow(2) or
# PhaseAbove(3):
#     if T >= ConstantValue  ->  phase = PhaseAbove (3)
#     else                   ->  phase = PhaseBelow (2)
# This mirrors that rule exactly using the (dimensional) marker temperature
# passed through the extended ABI.
#
# SCOPE: this only reproduces Check_Constant_Phase_Transition for
# number_phases = 1, PhaseDirection = BothWays, and no ResetParam (APS is
# left untouched here; the built-in Constant transition can also reset APS
# to 0 when ResetParam=APS is set - PT0 in the t16 .dat does not set it, so
# this is not exercised). It does NOT generalise to number_phases > 1,
# PhaseDirection = BelowToAbove/AboveToBelow, or ResetParam.
#
# Comparison note: because LaMEM's built-in Phase_Transition() applies its
# phase transitions PT0..PT3 in sequence within ONE call (and, in the full
# t16 .dat, PT2/PT3 act on phase 3, which PT0 itself can just have produced,
# i.e. those transitions cascade within the same time step), replacing only
# PT0 with this plugin while PT1-3 stay built-in would NOT be expected to
# give identical results to the fully-built-in run: the plugin call happens
# AFTER Phase_Transition() returns, so a marker that PT0 flips to phase 3
# this step would only become visible to the built-in PT2 (Clapeyron, 3<->5)
# on the FOLLOWING step, not immediately. The comparison actually run for
# this report therefore uses a copy of the .dat with ONLY PhaseTransition
# ID 0 present (PT1-3 removed), so built-in vs. plugin must match exactly
# with no ordering effects. See doc/phase_transition_plugin_PHASE1_REPORT.md.

const CONSTANT_VALUE = 1200.0   # dimensional, scal units
const PHASE_BELOW     = 2
const PHASE_ABOVE     = 3

Base.@ccallable function lamem_phase_transition(
        n::Csize_t,
        x::Ptr{Cdouble}, y::Ptr{Cdouble}, z::Ptr{Cdouble},
        T::Ptr{Cdouble}, p::Ptr{Cdouble}, time::Cdouble,
        sxx::Ptr{Cdouble}, syy::Ptr{Cdouble}, szz::Ptr{Cdouble},
        sxy::Ptr{Cdouble}, sxz::Ptr{Cdouble}, syz::Ptr{Cdouble},
        j2_stress_cell::Ptr{Cdouble}, j2_strainrate_cell::Ptr{Cdouble},
        eta_cell::Ptr{Cdouble}, aps_cell::Ptr{Cdouble},
        phase_in::Ptr{Cint}, phase_out::Ptr{Cint})::Cint
    try
        nn = Int(n)
        TT   = unsafe_wrap(Array, T, nn)
        Pin  = unsafe_wrap(Array, phase_in, nn)
        Pout = unsafe_wrap(Array, phase_out, nn)

        _ = (x, y, z, p, time, sxx, syy, szz, sxy, sxz, syz,
             j2_stress_cell, j2_strainrate_cell, eta_cell, aps_cell)

        changed = 0
        @inbounds for i in 1:nn
            ph = Pin[i]
            if ph == PHASE_BELOW || ph == PHASE_ABOVE
                newph = TT[i] >= CONSTANT_VALUE ? Cint(PHASE_ABOVE) : Cint(PHASE_BELOW)
                if newph != ph
                    changed += 1
                end
                ph = newph
            end
            Pout[i] = ph
        end
        return Cint(changed)
    catch
        # See ptlib.jl for why no I/O happens here under --trim=safe.
        return Cint(-1)
    end
end

end # module
