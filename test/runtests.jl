# This is the LaMEM testing framework, run through julia
#

using Test, GeophysicalModelGenerator

include("test_utils.jl")

@testset "LaMEM Testsuite" begin

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


# t3_SubductionMATLABinput - to be added  & changed to julia.
@testset "t3_Subduction" begin
    dir = "t3_SubductionMATLABinput";

    # input script 
    include(joinpath(dir,"CreateMarkers_Subduction.jl"));      

    ParamFile = "Subduction_MATLAB_Particles.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,), (rtol=1e-5,), (rtol=5e-4,));
    
    # test on 1 core
    # t3_Sub1_MATLAB_a_Direct_opt
    CreateMarkers_Subduction(dir, ParamFile, NumberCores=1)

    @test perform_lamem_test(dir,ParamFile,"Sub1_MATLAB_a_Direct_opt-p1.expected", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true)


    # t3_Sub1_MATLAB_b_MUMPS_opt                            
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=1e-5), (rtol=1e-5,atol=1e-5), (rtol=2.5e-4,atol=1e-3));
    
    ParamFile = "Subduction_MATLAB_Particles.dat";
    CreateMarkers_Subduction(dir, ParamFile, NumberCores=4)
    @test perform_lamem_test(dir,ParamFile,"Sub1_MATLAB_b_MUMPS_opt-p4.expected", 
                                keywords=keywords, accuracy=acc, cores=4, opt=true)
                        

    # t3_Sub1_MATLAB_c_MUMPS_deb                                 
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=2e-6), (rtol=1e-5,atol=3e-6), (rtol=2.5e-4,atol=3e-4));
    
    ParamFile = "Subduction_MATLAB_Particles4.dat";
    CreateMarkers_Subduction(dir, ParamFile, NumberCores=4)
    @test perform_lamem_test(dir,ParamFile,"Sub1_MATLAB_c_MUMPS_deb-p4.expected", 
                                args="-jp_pc_factor_mat_solver_type mumps",
                                keywords=keywords, accuracy=acc, cores=4, deb=true)
                        
    # t3_Sub1_MATLAB_d_MUMPS_MG_VEP_opt                                 
    # NOTE: This employs 1D grid refinement which does not work yet in julia
    #keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    #acc      = ((rtol=1e-6,atol=1e-6), (rtol=1e-5,atol=3e-6), (rtol=2.5e-4,atol=1e-4));
    
    #ParamFile = "Subduction_VEP.dat";
    #CreateMarkers_Subduction(dir, ParamFile, NumberCores=8)
    #@test perform_lamem_test(dir,ParamFile,"Sub1_MATLAB_d_MUMPS_MG_VEP_opt-p8.expected", 
    #                           args="",
    #                            keywords=keywords, accuracy=acc, cores=8, opt=true)
                        
end

@testset "t4_Localisation" begin
    dir = "t4_Loc";
    ParamFile = "localization.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-10), (rtol=1e-5,atol=2e-10), (rtol=1e-4,atol=1e-7));
    
    # Perform tests
    # t4_Loc1_a_MUMPS_VEP_opt
    @test perform_lamem_test(dir,"localization.dat","Loc1_a_MUMPS_VEP_opt-p4.expected",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=4, opt=true)

    # t4_Loc1_b_MUMPS_VEP_Reg_opt
    @test perform_lamem_test(dir,"localization_eta_min_reg.dat","Loc1_b_MUMPS_VEP_Reg_opt-p4.expected",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=4, opt=true)

    # t4_Loc1_c_Direct_VEP_opt                            
    @test perform_lamem_test(dir,"localization.dat","Loc1_c_Direct_VEP_opt-p1.expected",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true)
end

@testset "t5_Permeability" begin
    dir = "t5_Perm";
    ParamFile = "Permea.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,), (rtol=1e-5,), (rtol=1e-2,));
    
    # t5_Permeability_Direct_opt
    @test perform_lamem_test(dir,ParamFile,"Permeability_direct_opt-p4.expected", 
                            keywords=keywords, accuracy=acc, cores=4, opt=true)
end

@testset "t6_AdjointGradientScalingLaws_p2" begin
    dir = "t6_AdjointGradientScaling";
    
    keywords   = (  "|Div|_inf",
                    "|Div|_2",
                    "|mRes|_2",
                    "|                  rho[  0]",
                    "|                  eta[  1]",
                    "|                  eta[  0]",
                    "|   Prefactor A               :",
                    "|   Velocity check            :")

    acc        = (  (rtol=1e-7, atol=1e-6), 
                    (rtol=1e-5, atol=1e-5), 
                    (rtol=1e-4, atol=1e-5), 
                    (rtol=1e-6, atol=1e-6),
                    (rtol=1e-6, atol=1e-6),
                    (rtol=2e-8, atol=1e-6),
                    (rtol=2e-8, atol=1e-6),
                    (rtol=2e-8, atol=1e-6)
                );
    split_sign = ("=","=","=","","","","","")

    # Perform tests
    # t6_AdjointGradientScalingLaws_p2
    ParamFile = "t6_RTI_ScalingLaw.dat";
    @test perform_lamem_test(dir,ParamFile,"t6_AdjointGradientScaling_p2.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, split_sign=split_sign)

    # t6_AdjointGradientScalingLaws_SoftFilm
    ParamFile = "t6_RTI_ScalingLaw.dat";
    @test perform_lamem_test(dir,ParamFile,"t6_AdjointGradientScaling_SoftFilm_p1.expected",
                            args = "-surf_level 0.1 -eta[0] 10 -eta[1] 1 -coord_x -0.4,0.4 -FreeSurf_Wavelength 0.8", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, split_sign=split_sign)

end


@testset "t7_AdjointGradientInversion" begin
    dir = "t7_AdjointGradientInversion";
    
    # t7_AdjointGradientInversion_1
    keywords   = (  "| 1 Diff parameter value =",
                    "| 2 Diff parameter value =",
                    "| 1 Parameter value =",
                    "| 2 Parameter value =")

    acc        = (  (rtol=1e-6, atol=1e-6), 
                    (rtol=1e-5, atol=1e-5), 
                    (rtol=1e-5, atol=1e-5), 
                    (rtol=1e-5, atol=1e-6),
                    (rtol=1e-5, atol=1e-6),
                );

    # Perform tests
    ParamFile = "t7_Subduction2D_FreeSlip_Inversion.dat";
    @test perform_lamem_test(dir,ParamFile,"t7_AdjointGradientInversion_1.expected",
                            args="-nel_z 16 -nel_x 64",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # t7_AdjointGradientInversion_2
    keywords   = (  "| misfit          ",
                    "| misfit / misfit0",
                    "|   1 eta[0] =",
                    "|   2 rho[1] =",
                    "|  adjoint     1:   log10  eta[ 0]",
                    "|  adjoint     2:          rho[ 1]")

    acc        = (  (rtol=1e-5, atol=1e-6), 
                    (rtol=1e-5, atol=1e-5), 
                    (rtol=1e-5, atol=1e-5), 
                    (rtol=1e-5, atol=1e-6),
                    (rtol=1e-5, atol=1e-6),
                    (rtol=1e-5, atol=1e-6),
                );
    split_sign = ("=","=","","","","")

    ParamFile = "t7_Subduction2D_FreeSlip_Inversion.dat";
    @test perform_lamem_test(dir,ParamFile,"t7_AdjointGradientInversion_2.expected",
                            args="-tao_fmin 1e-6 -nel_z 16 -nel_x 64 -Inversion_EmployTAO 1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true); 

    # t7_AdjointGradientInversion_3
    keywords   = (  "| misfit          ",
                    "| misfit / misfit0",
                    "|   1 eta[0] =",
                    "|   2 rho[1] =",
                    "|       FD     1:   log10  eta[ 0]",
                    "|       FD     2:          rho[ 1]")

    acc        = (  (rtol=1e-3, atol=1e-6), 
                    (rtol=1e-3, atol=1e-5), 
                    (rtol=1e-3, atol=1e-5), 
                    (rtol=1e-5, atol=1e-6),
                    (rtol=1e-3, atol=1e-6),
                    (rtol=2e-4, atol=1e-6),
                );

    ParamFile = "t7_Subduction2D_FreeSlip_Inversion_FD.dat";
    @test perform_lamem_test(dir,ParamFile,"t7_AdjointGradientInversion_3.expected",
                            args="-tao_fmin 1e-6 -nel_z 16 -nel_x 32",
                            keywords=keywords, accuracy=acc, cores=1, opt=true) 

    # PSD paper inversion for nonlinear materials:
    # t7_PSDInversion_1
    keywords   = (  "| LS factor for 1.Parameter = ",
                    "|    F =",
                    "| 1 Parameter value = "
                 )

    acc        = (  (rtol=1e-1, atol=1e-6), 
                    (rtol=1e-3, atol=1e-5), 
                    (rtol=1e-5, atol=1e-5),
                );
    
    ParamFile = "t7_PSDInversionPaper.dat";
    @test perform_lamem_test(dir,ParamFile,"t7_PSDInversionPaper_1.expected",
                            args="-Inversion_rtol 4.6e-2",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)

    # PSD paper inversion for linear materials:   
    # t7_PSDInversion_2                            
    keywords   = (  "| LS factor for 1.Parameter = ",
                    "|    F =",
                    "| 1 Parameter value = "
                 )

    acc        = (  (rtol=1e-1, atol=1e-6), 
                    (rtol=1e-3, atol=1e-5), 
                    (rtol=5e-5, atol=1e-5),
                );

    ParamFile = "t7_PSDInversionPaper.dat";
    @test perform_lamem_test(dir,ParamFile,"t7_PSDInversionPaper_2.expected",
                            args="-nel_x 8 -nel_y 8 -nel_z 8  -n[0] 1 -n[1] 1 -n[2] 1  -Value[0] 135",
                            keywords=keywords, accuracy=acc, cores=2, opt=true) 

end

@testset "t8_AdjointGradients" begin
  include("AdjointGradients.jl")
end

@testset "t9_PhaseDiagrams" begin
    dir = "t9_PhaseDiagrams";
    ParamFile = "test_9_FallingBlock_PhaseDiagrams.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=1e-9));
    
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

@testset "t11_subgrid" begin
    dir = "t11_subgrid";
    ParamFile = "FallingBlock_mono_CoupledMG_RedundantCoarse.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-5,atol=1e-6), (rtol=1e-5, atol=1e-5), (rtol=1e-4,atol=1e-5));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"t11_Subgrid_opt-p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

end


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

# t13_Rheology0D/

# t14_1DStrengthEnvelope/

# t15_RTI/

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

    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-9));
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
    # t17_InflowOutflow2D_opt
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_2D.dat","InflowOutflow-2D_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)
    
    # 3D test
    # t17_InflowOutflow3D_opt
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-10), (rtol=1e-4,atol=2e-10));
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_3D.dat","InflowOutflow-3D_p4.expected",
                            keywords=keywords, accuracy=acc, cores=4, opt=true)

    # Test inflow/outflow conditions in 2D using optimized LaMEM   
    # t17_InflowOutflow2D_Pres_opt 
    acc      = ((rtol=2e-7,atol=1e-10), (rtol=1e-5, atol=1e-10), (rtol=1e-4,atol=2e-10), (rtol=1e-6,atol=1e-9));
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_2D_Perm.dat","InflowOutflow-2D_Perm_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # test_3D_Pres():
    # t17_InflowOutflow3D_Pres_opt
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11), (rtol=1e-7,atol=1e-11));
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_3D_Perm.dat","InflowOutflow-3D_Perm_p4.expected",
                            keywords=keywords, accuracy=acc, cores=4, opt=true)
                              
end

@testset "t18_SimpleShear" begin
    dir = "t18_SimpleShear";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-8), (rtol=1e-5, atol=2e-8), (rtol=1e-4,atol=2e-5));

    # test_xz
    @test perform_lamem_test(dir,"SS.dat","SimpleShear_xz-p2.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)
    
    # test_yz
    @test perform_lamem_test(dir,"SS.dat","SimpleShear_yz-p2.expected",
                            args="-exz_strain_rates 0 -eyz_strain_rates 1e-15 -eyz_num_periods 1",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)
    # test_xy
    @test perform_lamem_test(dir,"SS.dat","SimpleShear_xy-p2.expected",
                            args="-exz_strain_rates 0 -exy_strain_rates 1e-15 -exy_num_periods 1",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)

    # test_xz_yz
    @test perform_lamem_test(dir,"SS.dat","SimpleShear_xz_yz-p2.expected",
                            args="-exz_strain_rates 1e-15 -eyz_strain_rates 1e-15 -eyz_num_periods 1",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)
end

@testset "t19_CompensatedInflow" begin
    dir = "t19_CompensatedInflow";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_a
    # t19_CompensatedInflow
    @test perform_lamem_test(dir,"CompensatedInflow_test_2D.dat","CompensatedInflow-p1.expected",
                            args="-nstep_max 10",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # test_b
    # t19_CompensatedInflow3D
    @test perform_lamem_test(dir,"CompensatedInflow_test_3D.dat","CompensatedInflow3D-p2.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)

end

@testset "t20_FSSA" begin
    dir = "t20_FSSA";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-5,atol=1e-5), (rtol=1e-5, atol=1e-5), (rtol=1e-4,atol=1e-4));

    # t20_FSSA_1_opt
    @test perform_lamem_test(dir,"RTI_FSSA.dat","RTI_FSSA_1-p1.expected",
                            args="-nstep_max 20 -nel_x 50 -nel_z 100",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

end

@testset "t21_Passive_Tracer" begin
    dir = "t21_Passive_Tracer";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-9), (rtol=1e-5, atol=1e-8), (rtol=1e-4,atol=1e-4));

    # test_a
    # t21_Passive_Tracer_Always
    @test perform_lamem_test(dir,"Passive_tracer_ex2D.dat","Passive_tracer-2D_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # test_b
    # t21_Passive_Tracer_Condition
    @test perform_lamem_test(dir,"Passive_tracer_ex2D_Condition.dat","Passive_tracer-2D_Condition_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)
end

@testset "t22_RidgeGeom" begin
    dir = "t22_RidgeGeom";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_2D
    @test perform_lamem_test(dir,"ridge_geom_2D.dat","RidgeGeom2D.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # Test oblique ridge geometry conditions in 3D using optimized LaMEM   
    @test perform_lamem_test(dir,"ridge_geom_oblique_2cores.dat","RidgeGeom_oblique_2cores.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)
end

@testset "t23_Permeable" begin
    dir = "t23_Permeable";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=2e-9));

    # test_a
    @test perform_lamem_test(dir,"Permeable.dat","Permeable_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)
end

#=
MATLAB setup
@testset "t24_Erosion_Sedimentation" begin
    dir = "t24_Erosion_Sedimentation";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_a
    @test perform_lamem_test(dir,"Permeable.dat","Permeable_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

  
end
=#

@testset "t25_APS_Healing" begin
    dir = "t25_APS_Healing";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_2D
    @test perform_lamem_test(dir,"APS_Healing2D.dat","APS_Healing2D.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)
  
end


@testset "t26_Dike" begin
    dir = "t26_Dike";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_M1_2D
    @test perform_lamem_test(dir,"dike_M1_2D.dat","dike_M1_2D.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # test_M075_2D_2cores
    @test perform_lamem_test(dir,"dike_M075_2D_2cores.dat","dike_M075_2D_2cores.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)


    # test_variableM
    @test perform_lamem_test(dir,"dike_variableM.dat","dike_variableM.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # heat_kfac
    keywords = ("|eRes|_2",)
    acc      = ((rtol=1e-4,atol=1e-11),);
    @test perform_lamem_test(dir,"dike_heating_kfac.dat","dike_heating_kfac.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # heat_rhoA
    keywords = ("|eRes|_2",)
    acc      = ((rtol=1e-4,atol=6e-9),);
    @test perform_lamem_test(dir,"dike_heating_rhoA.dat","dike_heating_rhoA.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

end


@testset "t27_T-dep_Conductivity" begin
    dir = "t27_T-dep_Conductivity";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_2fields_dike():
    @test perform_lamem_test(dir,"t27_TDep_NuK_Conductivity.dat","t27_TDep_NuK_Conductivity.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

end

@testset "t28_HeatRecharge" begin
    dir = "t28_HeatRecharge";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=5e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-5,atol=1e-11));

    # test_recharge1
    @test perform_lamem_test(dir,"FallingBlockHeatReacharge1.dat","t28_HeatRecharge1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)

    # test_recharge2
    acc      = ((rtol=3e-6,atol=5e-7), (rtol=1e-5, atol=5e-6), (rtol=3e-5,atol=2e-5));
    @test perform_lamem_test(dir,"FallingBlockHeatReacharge2.dat","t28_HeatRecharge2.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true)
end


@testset "t29_PermeableSides_VelBoxes" begin
    dir = "t29_PermeableSides_VelBoxes";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=5e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-5,atol=1e-11));

    # test_permeableSides_VelBoxes
    @test perform_lamem_test(dir,"VelBoxes_Permeable_sides.dat","t29_PermeableSides_VelBoxes.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true)

end

@testset "t30_Timestep_Schedule" begin
    dir = "t30_Timestep_Schedule";
    keywords = ("Actual time step",)
    acc      = ((rtol=1e6,atol=1e-11),);

    # test_TS_Schedule():
    @test perform_lamem_test(dir,"TS_Schedule.dat","t30_TS_Schedule.expected",
                            keywords=keywords, accuracy=acc, cores=4, opt=true, split_sign=":")

end


end

