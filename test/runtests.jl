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


# t3_SubductionMATLABinput - to be added  

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