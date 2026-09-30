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
#
# BIT-EXACT REPRODUCTION VIA INTERNAL UNITS: Set_Constant_Phase_Transition
# (src/phase_transition.cpp) does NOT compare dimensional T against a
# dimensional ConstantValue at runtime. It non-dimensionalises ConstantValue
# ONCE, when the .dat file is read:
#
#     ph->ConstantValue = (ph->ConstantValue + scal->Tshift) / scal->temperature;
#
# and then compares the marker's raw INTERNAL P->T against that internal
# threshold every time step:
#
#     if (P->T >= PhaseTrans->ConstantValue) { ph = PH2; ... }
#
# This module reproduces that exact arithmetic path, not just an
# algebraically-equivalent dimensional comparison: CONSTANT_VALUE_INTERNAL
# below is computed with the SAME formula LaMEM itself applies to its
# ConstantValue (using the scaling struct's Tshift/temperature, read once
# per call since they are constant for a given run), and the comparison is
# performed against `m.T` after re-deriving the internal T via
# LaMEMPlugin.nondimensionalize_T - i.e. this plugin recomputes and compares
# in INTERNAL units, exactly mirroring LaMEM's own comparison, rather than
# comparing the wrapper's dimensional `m.T` against a dimensional 1200.
# (Algebraically the two are identical for finite floating-point values -
# dividing both sides of an inequality by the same positive `temperature`
# and adding/subtracting the same `Tshift` does not change its truth value
# except possibly at the exact boundary under floating-point rounding - but
# comparing in internal units, the way this file does, removes even that
# theoretical possibility of a rounding-induced mismatch, since it performs
# the identical sequence of floating-point operations LaMEM's own C code
# does.)
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
# on the FOLLOWING step, not immediately. The test therefore uses a copy of
# the .dat with ONLY PhaseTransition ID 0 present (PT1-3 removed), so built-in
# vs. plugin must match exactly with no ordering effects.

include("LaMEMPlugin.jl")
using .LaMEMPlugin

const CONSTANT_VALUE_DIM = 1200.0   # dimensional, scal units (Celsius in geo mode)
const PHASE_BELOW = Cint(2)
const PHASE_ABOVE = Cint(3)

function constant_transition_rule(m::MarkerView, constant_value_internal::Float64)
    if m.phase == PHASE_BELOW || m.phase == PHASE_ABOVE
        # compare in INTERNAL units, exactly like Check_Constant_Phase_Transition:
        # m.T_internal is LaMEM's raw P->T (unconverted), and
        # constant_value_internal was derived with the identical
        # (value + Tshift)/temperature formula LaMEM applies to its own
        # ConstantValue when the .dat is read - so this performs the exact
        # same floating-point comparison LaMEM's C code does.
        newph = m.T_internal >= constant_value_internal ? PHASE_ABOVE : PHASE_BELOW
        return (newph, m.T) # T itself is left untouched by this transition
    end
    return (m.phase, m.T)
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
        s = unsafe_load(scaling)
        constant_value_internal = (CONSTANT_VALUE_DIM + s.Tshift) / s.temperature

        rule = m -> constant_transition_rule(m, constant_value_internal)

        return lamem_pt_wrapper(rule, n, x, y, z, T, p, time,
            sxx, syy, szz, sxy, sxz, syz,
            j2_stress_cell, j2_strainrate_cell, eta_cell, aps_cell,
            phase_in, phase_out, T_out, scaling)
    catch
        return Cint(-1)
    end
end

Base.@ccallable function lamem_plugin_abi_version()::Cint
    return LaMEMPlugin.ABI_VERSION
end

end # module
