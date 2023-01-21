# This is the LaMEM testing framework, run through julia
#

using Test, GeophysicalModelGenerator

include("test_utils.jl")


@testset "t1_FB1_Direct" begin
    dir = "t1_FB1_Direct";
    ParamFile = "FallingBlock_mono_PenaltyDirect.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,), (rtol=1e-5,), (rtol=1e-4,));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"FB1_a_Direct_opt-p1.expected", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    @test perform_lamem_test(dir,ParamFile,"FB1_b_Direct_deb-p1.expected", 
                            keywords=keywords, accuracy=acc, cores=1, deb=true)

    @test perform_lamem_test(dir,ParamFile,"FB1_c_MUMPS_opt-p2.expected", 
                            keywords=keywords, accuracy=acc, cores=2, opt=true,
                            args="-jp_pc_factor_mat_solver_package mumps")

    #@test perform_lamem_test(dir,ParamFile,"FB1_d_PaStiX_opt-p4.expected", 
    #                        keywords=keywords, accuracy=acc, cores=4, opt=true,
    #                        args="-jp_pc_factor_mat_solver_package pastix")

end

@testset "t2_FB2_MG" begin
    dir = "t2_FB2_MG";
    ParamFile = "FallingBlock_mono_CoupledMG_RedundantCoarse.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,), (rtol=1e-5,), (rtol=1e-4,));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"FB2_a_CoupledMG_opt-p1.expected", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

end


# t3_SubductionMATLABinput - to be added  & changed to julia


@testset "t4_Loc" begin
    dir = "t4_Loc";
    ParamFile = "localization.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-10), (rtol=1e-5,atol=1e-10), (rtol=1e-4,));
    
    # Perform tests
    @test perform_lamem_test(dir,"localization.dat","Loc1_a_MUMPS_VEP_opt-p4.expected",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=4, opt=true)

    @test perform_lamem_test(dir,"localization_eta_min_reg.dat","Loc1_b_MUMPS_VEP_Reg_opt-p4.expected",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=4, opt=true)
end

@testset "t5_Perm" begin
    dir = "t5_Perm";
    ParamFile = "Permea.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"Permeability_direct_opt-p4.expected",
                            keywords=keywords, accuracy=acc, cores=4, opt=true)

end


#=
# This does not work currently
@testset "t6_AdjointGradientScaling" begin
    dir = "t6_AdjointGradientScaling";
    ParamFile = "t6_RTI_ScalingLaw.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2", "|   Prefactor A               :")
    acc      = ((rtol=1e-7,), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11), (rtol=2e-8,));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"t6_AdjointGradientScaling_p2.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)

    @test perform_lamem_test(dir,ParamFile,"t6_AdjointGradientScaling_SoftFilm_p1.expected",
                            args="-surf_level 0.1 -eta[0] 10 -eta[1] 1 -coord_x -0.4,0.4 -FreeSurf_Wavelength 0.8",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)


end
=#

#=
@testset "t7_AdjointGradientInversion" begin

    dir = "t7_AdjointGradientInversion";
    ParamFile = "t7_Subduction2D_FreeSlip_Inversion.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2", "|   Prefactor A               :")
    acc      = ((rtol=1e-7,), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11), (rtol=2e-8,));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"t6_AdjointGradientScaling_p2.expected",
                            args="-nel_z 16 -nel_x 64",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)


end
=#

#t8 to be done

@testset "t9_PhaseDiagrams" begin
    dir = "t9_PhaseDiagrams";
    ParamFile = "test_9_FallingBlock_PhaseDiagrams.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=1e-11));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"test_9_FallingBlock_PhaseDiagrams.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)

end

#=
# this ia a more complicated one, that requires a devoted script (with plotting)

@testset "t10_Compressibility" begin
    dir = "t10_Compressibility";
    ParamFile = "test_9_FallingBlock_PhaseDiagrams.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=1e-11), (rtol=2e-8,));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"test_9_FallingBlock_PhaseDiagrams.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)

end
=#

#=
fails
@testset "t11_subgrid" begin
    dir = "t11_subgrid";
    ParamFile = "FallingBlock_mono_CoupledMG_RedundantCoarse.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=1e-11));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"t11_Subgrid_opt-p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

end
=#

#=
another more complicated test
@testset "t12_Temperature_diffusion" begin
    dir = "t12_Temperature_diffusion";
    ParamFile = "FallingBlock_mono_CoupledMG_RedundantCoarse.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=1e-11));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"t11_Subgrid_opt-p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

end
=#


@testset "t16_PhaseTransitions" begin
    dir = "t16_PhaseTransitions";
    ParamFile = "Plume_PhaseTransitions.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=1e-9));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"PhaseTransitions-p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    @test perform_lamem_test(dir,ParamFile,"PhaseTransitions-FreeSlip_p1.expected",
                            args="-open_top_bound 0 -act_press_shift 1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    @test perform_lamem_test(dir,"Plume_PhaseTransitions_Melting.dat","PhaseTransitions-Melting_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # Tests phase transitions with X/Z and Box coordinates
    @test perform_lamem_test(dir,"Plume_PhaseTransitions_Box_XZ.dat","PhaseTransitions-XBox-p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # Tests phase transition triggered by time       
    @test perform_lamem_test(dir,"TimeTransition.dat","TimeTransition-p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)
                            
    # Test dike feature using optimized LaMEM
    @test perform_lamem_test(dir,"PhaseTransitionBox_move.dat","PhaseTransitionBox_move.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

end


@testset "t17_InflowOutflow" begin
    dir = "t17_InflowOutflow";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2","|eRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11), (rtol=1e-7,atol=1e-11));
    
    # 2D test
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_2D.dat","InflowOutflow-2D_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)
    # 3D test
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_3D.dat","InflowOutflow-2D_p1.expected",
                            keywords=keywords, accuracy=acc, cores=4, opt=true)

    # Test inflow/outflow conditions in 2D using optimized LaMEM    

end