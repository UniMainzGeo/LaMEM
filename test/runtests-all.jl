# This is the LaMEM testing framework, run through julia
#
# Start it from 
using LaMEM_C
using Test
using GeophysicalModelGenerator
using PETSc_jll

const create_plots = false
if create_plots
    using CairoMakie
end

# Read all julia IO functions
include("julia/IO_functions.jl")  # copied from LaMEM.jl; we do not want to make LaMEM.jl a depencency here as it fixes the PETSc_jll version
using .IO_functions

include("julia/run_lamem_save_grid.jl")  

if "use_dynamic_lib" in ARGS
    global use_dynamic_lib=true
else
    global use_dynamic_lib=false
end
  #global use_dynamic_lib=true
test_mumps=true # if we do this later on windows, we have to deactivate this

if "is64bit" in ARGS
    global is64bit=true
else
    global is64bit=false
end

@show use_dynamic_lib create_plots
include("test_utils.jl")        # test-framework specific functions

test_dir = pwd()

#---------------------------------------------------------------------------
maintenance = false # set to true when designing/debugging tests

if maintenance

	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	# WARNING! HANDLE WITH CARE
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	update_expected = false # set to true when ready to update
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	if update_expected
	
		print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
		print("WARNING! YOU ARE ABOUT TO OVERWRITE THE EXPECTED FILES\n")
		print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
	
		clean_files = true # clean files when updating expected
	else
		clean_files = false # no noeed to clean files, keep working
	end
	
else
	# normal mode (do not update the expected files, just clean the output/work files)
	update_expected = false 
	clean_files     = true
end
#---------------------------------------------------------------------------
@testset "LaMEM Testsuite" verbose=true begin
#---------------------------------------------------------------------------
@testset "t01_FB1_Direct" verbose=true begin
    cd(test_dir)
    dir = "t01_FB1_Direct";
    
    ParamFile = "FallingBlock_Direct_Default.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,), (rtol=1e-5,), (rtol=1e-4,));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"FB1_a_Direct_opt", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    @test perform_lamem_test(dir,ParamFile,"FB1_b_Direct_deb", 
                            keywords=keywords, accuracy=acc, cores=1, deb=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    ParamFile = "FallingBlock_Direct_MUMPS.dat";

    @test perform_lamem_test(dir,ParamFile,"FB1_c_MUMPS_opt", 
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t02_FB2_MG" begin

    cd(test_dir)
    dir = "t02_FB2_MG";
    
    ParamFile = "FallingBlock_mono_CoupledMG_RedundantCoarse.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,), (rtol=1e-5,), (rtol=1e-4,));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"FB2_a_CoupledMG_opt", 
                            keywords=keywords, accuracy=acc, cores=4, deb=true, opt=false, mpiexec=mpiexec, debug=false,
                                create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t03_Subduction" begin
    cd(test_dir)
    dir = "t03_SubductionGMGinput";
    
    # input script 
    include(joinpath(dir,"CreateMarkers_Subduction.jl"));      

    ParamFile = "Subduction_GMG_Particles_Default.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=5e-7), (rtol=1e-5,atol=1e-5), (rtol=5e-4,atol=1e-3));
    
    # test on 1 core
    # Sub1_a_Direct_opt
    CreateMarkers_Subduction(dir, ParamFile, NumberCores=1, mpiexec=mpiexec)

    @test perform_lamem_test(dir,ParamFile,"Sub1_a_Direct_opt", 
                            args="-nstep_max 2",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # Sub1_b_MUMPS_opt                            
    ParamFile = "Subduction_GMG_Particles_MUMPS.dat";

    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=1e-5), (rtol=1e-5,atol=1e-5), (rtol=2.5e-4,atol=1e-3));
        
    CreateMarkers_Subduction(dir, ParamFile, NumberCores=4, mpiexec=mpiexec, is64bit=is64bit)
    @test perform_lamem_test(dir,ParamFile,"Sub1_b_MUMPS_opt", 
                                args="-nstep_max 2",
                                keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                                create_expected_file=update_expected, clean_dir=clean_files)

    # Sub1_c_MUMPS_deb    
    ParamFile = "Subduction_GMG_Particles_MUMPS.dat";
                     
    # writing parallel marker files doesn't work in CI with 64 bit atm 
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=2e-6), (rtol=1e-5,atol=3e-6), (rtol=2.5e-4,atol=3e-4));
        
    CreateMarkers_Subduction(dir, ParamFile, NumberCores=4, mpiexec=mpiexec, is64bit=is64bit)
    @test perform_lamem_test(dir,ParamFile,"Sub1_c_MUMPS_deb", 
                                args="-nstep_max 2",
                                keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                                create_expected_file=update_expected, clean_dir=clean_files)
                        
    # Sub1_d_MUMPS_MG_VEP_opt                                 
    ParamFile = "Subduction_VEP.dat";
    # NOTE: This employs 1D grid refinement
    include(joinpath(dir,"CreateMarkers_SubductionVEP_parallel.jl"));      
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=1e-6), (rtol=1e-5,atol=3e-6), (rtol=2.5e-4,atol=1e-4));
        
    CreateMarkers_SubductionVEP(dir, ParamFile, NumberCores=2, mpiexec=mpiexec,  is64bit=is64bit)
    @test perform_lamem_test(dir,ParamFile,"Sub1_d_MUMPS_MG_VEP_opt", 
                                args="-nstep_max 2",
                                keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                                create_expected_file=update_expected, clean_dir=clean_files)    
end
#---------------------------------------------------------------------------
@testset "t04_Localisation" begin
    cd(test_dir)
    dir = "t04_Loc";
    
    ParamFile = "localization.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-10), (rtol=1e-5,atol=2e-9), (rtol=1e-4,atol=1e-7));
    
    # Perform tests
    if test_mumps & !is64bit
        # This test has issues on github actions with 64bit but works fine on our machines and with 32bit.

        # Loc1_a_MUMPS_VEP_opt
        @test perform_lamem_test(dir,"localization.dat","Loc1_a_MUMPS_VEP_opt",
                                args="-nstep_max 20", 
                                keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            	create_expected_file=update_expected, clean_dir=clean_files)

	end
    # Loc1_b_MUMPS_VEP_Reg_opt
    @test perform_lamem_test(dir,"localization_eta_min_reg.dat","Loc1_b_MUMPS_VEP_Reg_opt",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # Loc1_c_Direct_VEP_opt                            
    @test perform_lamem_test(dir,"localization.dat","Loc1_c_Direct_VEP_opt",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)


    # Loc1_d_MUMPS_VEP_VPReg_opt
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=5e-2,atol=1e-8), (rtol=2e-3,atol=5e-9), (rtol=2e-3,atol=2e-7));
    
    @test perform_lamem_test(dir,"localization_eta_vp_reg.dat","Loc1_d_MUMPS_VEP_VPReg_opt",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t05_Permeability" begin
    cd(test_dir)
    dir = "t05_Perm";
    
    ParamFile = "Permea.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,), (rtol=1e-5,), (rtol=1e-2,atol=1e-8));
    
    # Permeability_Direct_opt
    @test perform_lamem_test(dir,ParamFile,"Permeability_direct_opt", 
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t06_AdjointGradientScalingLaws_p2" begin
    cd(test_dir)
    dir = "t06_AdjointGradientScaling";
    
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
	# AdjointGradientScalingLaws_p2
    ParamFile = "RTI_ScalingLaw.dat";
    @test perform_lamem_test(dir,ParamFile,"AdjointGradientScaling_p2",
                             keywords=keywords, accuracy=acc, cores=2, opt=true, split_sign=split_sign, mpiexec=mpiexec,
                             create_expected_file=update_expected, clean_dir=clean_files)
        
    # AdjointGradientScalingLaws_SoftFilm
    ParamFile = "RTI_ScalingLaw.dat";
    @test perform_lamem_test(dir,ParamFile,"AdjointGradientScaling_SoftFilm_p1",
                            args = "-surf_level 0.1 -eta[0] 10 -eta[1] 1 -coord_x -0.4,0.4 -FreeSurf_Wavelength 0.8", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, split_sign=split_sign,  mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

end
#---------------------------------------------------------------------------
@testset "t07_AdjointGradientInversion" begin
    cd(test_dir)
    dir = "t07_AdjointGradientInversion";
    
    # AdjointGradientInversion_1
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
    ParamFile = "Subduction2D_FreeSlip_Inversion.dat";
    @test perform_lamem_test(dir,ParamFile,"AdjointGradientInversion_1",
                            args="-nel_z 16 -nel_x 64",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # AdjointGradientInversion_2
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

    ParamFile = "Subduction2D_FreeSlip_Inversion.dat";
    @test perform_lamem_test(dir,ParamFile,"AdjointGradientInversion_2",
                            args="-tao_fmin 1e-6 -nel_z 16 -nel_x 64 -Inversion_EmployTAO 1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files) 

    # AdjointGradientInversion_3
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

    ParamFile = "Subduction2D_FreeSlip_Inversion_FD.dat";
    @test perform_lamem_test(dir,ParamFile,"AdjointGradientInversion_3",
                            args="-tao_fmin 1e-6 -nel_z 16 -nel_x 32",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files) 

    # PSD paper inversion for nonlinear materials:
    # PSDInversion_1
    keywords   = (  "| LS factor for 1.Parameter = ",
                    "|    F =",
                    "| 1 Parameter value = "
                 )

    acc        = (  (rtol=1e-1, atol=1e-6), 
                    (rtol=1e-3, atol=1e-5), 
                    (rtol=1e-5, atol=1e-5),
                );
    
    ParamFile = "PSDInversionPaper.dat";
    @test perform_lamem_test(dir,ParamFile,"PSDInversionPaper_1",
                            args="-Inversion_rtol 4.6e-2 -nel_x -nel_y 8 -nel_z 8",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # PSD paper inversion for linear materials:   
    # PSDInversion_2                            
    keywords   = (  "| LS factor for 1.Parameter = ",
                    "|    F =",
                    "| 1 Parameter value = "
                 )

    acc        = (  (rtol=1e-1, atol=1e-6), 
                    (rtol=1e-3, atol=1e-5), 
                    (rtol=5e-5, atol=1e-5),
                );

    ParamFile = "PSDInversionPaper.dat";
    @test perform_lamem_test(dir,ParamFile,"PSDInversionPaper_2",
                            args="-nel_x 8 -nel_y 8 -nel_z 8  -n[0] 1 -n[1] 1 -n[2] 1  -Value[0] 135",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files) 
end
#---------------------------------------------------------------------------
@testset "t08_AdjointGradients" begin
	cd(test_dir)
	dir = "t08_AdjointGradients";
	# This tests various adjoint gradient cases
    
	# AdjointGradients_Sphere_ND_all
	keywords   = (  "|Div|_inf",
	               "|Div|_2",
	               "|mRes|_2",
	               "|           delta(rho)[  1]",
	               "|                  eta[  0]",
	               "|   Velocity check            :",
	               "|  adjoint     2:          eta[ 0]")
	
	acc        = (  (rtol=1e-7, atol=1e-6), 
	               (rtol=1e-8, atol=1e-5), 
	               (rtol=1e-8, atol=1e-5), 
	               (rtol=1e-6, atol=1e-5), 
	               (rtol=1e-6, atol=1e-5), 
	               (rtol=1e-6, atol=1e-5), 
	               (rtol=1e-6, atol=1e-5), 
	            );
	
	ParamFile = "AdjointGradients.dat";
	@test perform_lamem_test(dir,ParamFile,"AdjointGradients_Sphere_ND_all",
	                        args="",
	                        keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
	                        create_expected_file=update_expected, clean_dir=clean_files)
	
	# AdjointGradients_CompareGradients_1
	keywords   = (  "|Div|_inf",
	               "|Div|_2",
	               "|mRes|_2",
	               "|       FD     1:          eta[ 1]",
	               "|  adjoint     2:          eta[ 1]",
	               "|       FD     3:          eta[ 0]",
	               "|  adjoint     4:          eta[ 0]")
	
	acc        = (  (rtol=1e-7, atol=1e-6), 
	               (rtol=1e-8, atol=1e-5), 
	               (rtol=1e-8, atol=1e-5), 
	               (rtol=1e-8, atol=1e-5), 
	               (rtol=1e-8, atol=1e-5), 
	               (rtol=1e-8, atol=1e-5), 
	               (rtol=1e-8, atol=1e-5), 
	            );
	
	
	ParamFile = "AdjointGradients_CompareGradients.dat";
	@test perform_lamem_test(dir,ParamFile,"AdjointGradients_CompareGradients_1",
	                        args="",
	                        keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
	                        create_expected_file=update_expected, clean_dir=clean_files)
	
	# AdjointGradients_CompareGradients_geo
	keywords   = (  "|       FD     1:          eta[ 1]",
	                "|  adjoint     2:          eta[ 1]",
	                "|       FD     3:          eta[ 0]",
	                "|  adjoint     4:          eta[ 0]",
	                "|           delta(rho)[  1]")
	
	acc        = (  (atol=1e-30, ), 
	                (atol=1e-30, ), 
	                (atol=1e-28, ), 
	                (atol=1e-28, ), 
	                (atol=1e-3,  ),  
	             );
	
	ParamFile = "AdjointGradients_CompareGradients_geo.dat";
	@test perform_lamem_test(dir,ParamFile,"AdjointGradients_CompareGradients_geo",
	                        args="",
	                        keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
	                        create_expected_file=update_expected, clean_dir=clean_files)
	
	#AdjointGradients_CompareGradients_2
	keywords   = (  "|Div|_inf",
	                "|Div|_2",
	                "|mRes|_2",
	                "|       FD     1:            n[ 0]",
	                "|  adjoint     2:            n[ 0]",
	                "|   Prefactor A               :")
	
	acc        = (  (rtol=1e-7, atol=1e-6), 
	                (rtol=1e-8, atol=1e-5), 
	                (rtol=1e-8, atol=1e-5), 
	                (rtol=1e-8, atol=1e-5), 
	                (rtol=2e-6, atol=1e-5), 
	                (rtol=1e-8, atol=1e-5), 
	             );
	
	ParamFile = "AdjointGradients_CompareGradients_2.dat";
	@test perform_lamem_test(dir,ParamFile,"AdjointGradients_CompareGradients_2",
	                        args="",
	                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
	                        create_expected_file=update_expected, clean_dir=clean_files)
	
	# Adjoint_Subduction2D_FreeSlip
	keywords   = (  "|Div|_inf",
	                "|Div|_2",
	                "|mRes|_2",
	                "|                  eta[  0]",
	                "|           delta(rho)[  1]",
	                "|           delta(rho)[  2]",
	                "|                  eta[  2]",
	                "|                  eta[  1]",
	                "|      log10       eta[  0]",
	                "|       FD     7:           fr[ 2]"
	                )
	
	acc        = (  (rtol=1e-7, atol=1e-6), 
	                (rtol=1e-8, atol=1e-5), 
	                (rtol=1e-8, atol=1e-5), 
	                (rtol=1e-8, atol=1e-5), 
	                (rtol=1e-8, atol=1e-5), 
	                (rtol=1e-8, atol=1e-5),
	                (rtol=1e-8, atol=1e-5),
	                (rtol=1e-8, atol=1e-5),
	                (rtol=1e-3, atol=1e-5),
	                (rtol=1e-8, atol=1e-5),
	             );
	
	ParamFile = "Subduction2D_FreeSlip_DirectSolver.dat";
	@test perform_lamem_test(dir,ParamFile,"Subduction2D_FreeSlip_DirectSolver_p1",
	                        args="-nel_y 2",
	                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
	                        create_expected_file=update_expected, clean_dir=clean_files)
	
	# Adjoint_PSD
	keywords   = (  "|  adjoint     1:          rho[ 2]",
	                "|       FD     2:          rho[ 2]",
	                "| Current Cost function = "
	                )
	
	acc        = (  (rtol=1e-6, atol=1e-6), 
	                (rtol=1e-6, atol=1e-5), 
	                (rtol=1e-6, atol=1e-5), 
	             );
	
	ParamFile = "FB_PSDTest.dat";
	@test perform_lamem_test(dir,ParamFile,"FB_PSDTest_p1",
	                        args="-nel_x 8 -nel_y 8 -nel_z 8 ",
	                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
	                        create_expected_file=update_expected, clean_dir=clean_files)
	
	# Adjoint_rho_SensitivityKernel_PSD

	keywords   = ( "|Div|_inf",
	               "|Div|_2",
	               "|mRes|_2",
	               "| Current Cost function = "
	               )
	
	acc        = (  (rtol=1e-7, atol=1e-6), 
	               (rtol=1e-5, atol=1e-5), 
	               (rtol=1e-4, atol=1e-5), 
	               (rtol=1e-6, atol=1e-5), 
	            );   
	
	ParamFile = "AdjointGradients_SensitivityKernel_PSD.dat";
	@test perform_lamem_test(dir,ParamFile,"Adjoint_rho_SensitivityKernel_PSD_p2",
	                        args="",
	                        keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
	                        create_expected_file=update_expected, clean_dir=clean_files)

	
	# Adjoint_n_SensitivityKernel_PSD
	keywords   = ( "|   Norm of field gradient vector :",
	                )
	
	acc        = (  (rtol=5e-1, atol=1e-6), 
	             );
	split_sign       = (":",)       
	
	ParamFile = "PSDKernelPaper.dat";
	@test perform_lamem_test(dir,ParamFile,"Adjoint_n_SensitivityKernelPaper_PSD",
	                        args="-nel_x 8  -nel_y 8 -nel_z 8 ",
	                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
	                        create_expected_file=update_expected, clean_dir=clean_files)
	
	# Adjoint_eta0_SensitivityKernel_PSD
	keywords   = ( "|   Norm of field gradient vector :",
	                )
	
	acc        = (  (rtol=5e-1, atol=1e-6), 
	             );
	split_sign       = (":",)       
	
	ParamFile = "PSDKernelPaper.dat";
	@test perform_lamem_test(dir,ParamFile,"Adjoint_eta0_SensitivityKernelPaper_PSD",
	                        args="-nel_x 8  -nel_y 8 -nel_z 8 -Type[0] eta0",
	                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec, 
	                        split_sign=split_sign,
	                        create_expected_file=update_expected, clean_dir=clean_files)

end
#---------------------------------------------------------------------------
@testset "t09_PhaseDiagrams" begin
    cd(test_dir)
    dir = "t09_PhaseDiagrams";
    
    ParamFile = "FallingBlock_PhaseDiagrams.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-5,atol=1e-8));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"FallingBlock_PhaseDiagrams",
                            args="-mfmax 0.15",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t10_Compressibility" begin
	
	# this ia a more complicated one, that requires a devoted script (with plotting)

    cd(test_dir)
    dir = "t10_Compressibility";
    
    ParamFile = "Compressible1D_withSaltandBasement.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-7), (rtol=1e-5, atol=1e-7), (rtol=2e-6,atol=1e-4), (rtol=2e-8,atol=1e-6));
    
    # Perform tests

    # test_a -------
    @test perform_lamem_test(dir,ParamFile,"Compressibility_opt",
                            keywords=keywords, accuracy=acc, cores=1, opt=true,
                            create_expected_file=update_expected, clean_dir=false)
    
    # load the data
    data, t = read_LaMEM_timestep("output", 20, dir, last=true);

    # extract 1D profiles
    include(joinpath(dir,"Compressibility_analytics.jl"))
    phase_vec,ρ, z, Szz_vec, Sxx_vec, Pf_vec, τII_vec = extract_1D_profiles(data, dir)

    # 1D analytical solution
    Sv_a, Pf_a, P_hydro_a, Sh_a = AnalyticalSolution(ρ, phase_vec, z)

    # Compute difference with analytical solution
    @test norm(Szz_vec - Sv_a) ≈ 1.0769790188863786 rtol=1e-3
    @test norm(Sxx_vec - Sh_a) ≈ 19.596753495502448 rtol=1e-4
    @test norm(Pf_vec - Pf_a) ≈ 4.676818965337232 rtol=1e-5
    
    # Create plot with stress & analytical solution
    #Plot_vs_analyticalSolution(data, dir,"Compressible1D_output_1Core.png")
    if clean_files
    	clean_directory(dir)
    end
    # --------------

end
#---------------------------------------------------------------------------
@testset "t11_Subgrid" begin
    cd(test_dir)
    dir = "t11_Subgrid";
    
    ParamFile = "FallingBlockCoupledMG.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-1,atol=1e-5), (rtol=1e-5, atol=1e-5), (rtol=1e-4,atol=1e-5));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"Subgrid_opt",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                        	create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t12_Temperature_diffusion" begin
    cd(test_dir)
    dir = "t12_Temperature_diffusion";
    include(joinpath(dir,"Temp_setup.jl"))
    ParamFile = "Temperature_diffusion.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=5e-11));

    # ---
    # Perform tests
    CreateMarkers_Temperature(dir, "Temperature_diffusion.dat", "./markers_pT1"; NumberCores=1)

    @test perform_lamem_test(dir,ParamFile,"TpD_a",
                            args="-mark_load_file ./markers_pT1/mdb",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)
    
    data, t1 = read_LaMEM_timestep("t13", 1, dir); T1=data.fields.temperature[1,1,:]; 
    data, t3 = read_LaMEM_timestep("t13", 3, dir); T3=data.fields.temperature[1,1,:];
    data, t5 = read_LaMEM_timestep("t13", 5, dir); T5=data.fields.temperature[1,1,:];
    z = data.z.val[1,1,:]

    T_a5 = Analytical_1D(z, t5)
    @test norm(T_a5 - T5)/length(T5) ≈ 0.03356725901141689

    if create_plots
        Plot_Analytics_vs_Numerics(z,T_a5, T5, dir, "T_anal3.png")
    end
    if clean_files
    	clean_directory(dir)
    end
    # ---
   
    # halfspace cooling test ----
    ParamFile = "Temperature_diffusion_1Dhalfspace.dat"
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2","|T|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=5e-11),  (rtol=1e-1,atol=5e-11));

    @test perform_lamem_test(dir,ParamFile,"Temperature_diffusion",
                args="-printNorms 1",
                keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                create_expected_file=update_expected, clean_dir=clean_files)

	if clean_files
		clean_test_directory(dir)
	end 
    
    # ---
end
#---------------------------------------------------------------------------
# t13_Rheology0D
@testset "t13_Rheology0D" begin
    cd(test_dir)
    dir = "t13_Rheology0D";
    include(joinpath(dir,"Rheology0D.jl"))
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-9), (rtol=1e-6, atol=1e-9), (rtol=2e-6,atol=1e-9));

    # ---
    # Viscoelastic rheology
    @test perform_lamem_test(dir,"Rheology_VE_0D.dat","Rheology_VE_0D",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)

    # compare with analytics    
    FileName = "Rheolog0D_VE"                        
    t_vec, τII_LaMEM = StressTime_0D(FileName, dir);
    τII_anal = Viscoelastoplastic0D(5e10, 1e22, 1e-15, t_vec);      
    @test norm(τII_LaMEM-τII_anal/1e6)/length(τII_LaMEM) ≈ 0.12480014617816898  rtol = 1e-4

    # Create plot
    if create_plots
        t_anal = range(0,t_vec[end],200)
        τII_anal1 = Viscoelastoplastic0D(5e10, 1e22, 1e-15, t_anal)
        Plot_StressStrain(t_anal,τII_anal1/1e6, t_vec, τII_LaMEM, dir, "Viscoelastic0D.png")
        
        if clean_files
    		clean_directory(dir)
    	end
    end
    # ---

    # ---
    # Viscoelastoplastic rheology
    @test perform_lamem_test(dir,"Rheology_VEP_0D.dat","Rheology_VEP_0D",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)

    # compare with analytics     
    FileName = "Rheolog0D_VEP"
    t_vec, τII_LaMEM = StressTime_0D(FileName, dir);
    YieldStress = 1e7  
    τII_anal = Viscoelastoplastic0D(5e10, 1e22, 1e-15, t_vec, YieldStress);    
    @test norm(τII_LaMEM-τII_anal/1e6)/length(τII_LaMEM) ≈ 0.05341838341184021 rtol = 1e-4

    # Create plot
    if create_plots
        t_anal = range(0,t_vec[end],200)
        τII_anal1 = Viscoelastoplastic0D(5e10, 1e22, 1e-15, t_anal, YieldStress)
        Plot_StressStrain(t_anal,τII_anal1/1e6, t_vec, τII_LaMEM, dir, "Viscoelastoplastic0D.png")

        if clean_files
    		clean_directory(dir)
    	end
    end
    # ---

    # ---
    # Viscoelastoplastic rheology with nonlinear dislocation creep viscosity
    @test perform_lamem_test(dir,"Rheology_DislocationCreep_VE_0D.dat","Rheology_DislocationCreep_VE_0D",
                        keywords=keywords, accuracy=acc, cores=1, opt=true,  mpiexec=mpiexec,
                        create_expected_file=update_expected, clean_dir=false)
                        
    FileName = "Rheolog0D_DislocationCreep_VE"
    t_vec, τII_LaMEM = StressTime_0D(FileName, dir);
    YieldStress = 1e10  
    data,t = read_LaMEM_timestep(FileName, 0, dir, fields=("temperature [C]",));
    T = mean(data.fields.temperature)

    # Create plot
    if create_plots
        ε = 1e-15;
        t_anal, τII_anal1, τII_no_iter = Viscoelastoplastic0D_dislocationcreep(T, ε, maximum(t_vec))
        Plot_StressStrain(t_anal,τII_anal1/1e6, t_vec, τII_LaMEM, dir, "Viscoelastic0D_dislocationCreep.png", τII_no_iter=τII_no_iter/1e6)
    	
    	if clean_files
    		clean_directory(dir)
   		end
        
    end
    # ---
    
    # ---
    # Viscoelastoplastic rheology with nonlinear dislocation creep viscosity
    @test perform_lamem_test(dir,"Rheology_DislocationCreep_VEP_0D.dat","Rheology_DislocationCreep_VEP_0D",
                        keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                        create_expected_file=update_expected, clean_dir=false)
                        
    FileName = "Rheolog0D_DislocationCreep_VEP"
    t_vec, τII_LaMEM = StressTime_0D(FileName, dir);
    YieldStress = 15e6  
    data,t = read_LaMEM_timestep(FileName, 0, dir, fields=("temperature [C]",));
    T = mean(data.fields.temperature)

    # Create plot
    if create_plots
        t_anal, τII_anal1, τII_no_iter = Viscoelastoplastic0D_dislocationcreep(T, ε, maximum(t_vec), YieldStress)
        Plot_StressStrain(t_anal,τII_anal1/1e6, t_vec, τII_LaMEM, dir, "Viscoelastoplastic0D_dislocationCreep.png", τII_no_iter=τII_no_iter/1e6)
		
		if clean_files
			clean_directory(dir)
    	end

    end
    # ---

    # ---
    # Stress-strainrate for linear viscous rheologies
    ε = [-1e-13 -1e-14 -1e-15 -1e-16 -1e-17]
    FileName = "Rheology_linearViscous_0D.dat"
    τ = StressStrainrate0D_LaMEM(FileName, dir, "Rheolog0D_linearViscous", ε)
    slope = (log10.(-ε[end])-log10.(-ε[1]) )/(log10.(τ[end])-log10.(τ[1]))
    @test slope ≈ 1.0 rtol = 1e-6
    
    if create_plots
        τ_anal = -2*ε[:]*1e21/1e6
        Plot_StressStrainrate(ε, τ, τ_anal,  dir, "Stress_Strainrate_linearViscous.png")

		if clean_files
			clean_directory(dir)
    	end

    end
    # ---

    # ---
    # Stress-strainrate for dislocation creep rheologies
    FileName = "Rheology_PowerlawCreep_DryOlivine_0D.dat"
    τ = StressStrainrate0D_LaMEM(FileName, dir, "Rheolog0D_DryOlivine", ε)
    slope = (log10.(-ε[end])-log10.(-ε[1]) )/(log10.(τ[end])-log10.(τ[1]))
    @test slope ≈ 3.5 rtol=1e-2 

    # add analytical solution for DC
   
    T=1000;
    τ_anal = AnalyticalSolution_DislocationCreep("DryOlivine", T, ε)/1e6
    @test norm(τ_anal[:] .- τ[:]) ≈ 0.2009862117696578 rtol = 1e-4

    if create_plots
        Plot_StressStrainrate(ε, τ, τ_anal,  dir, "Stress_Strainrate_DryOlivine_DC.png")
        
        # clear all files in the test directory
        
		if clean_files
			clean_test_directory(dir)
		end         
    end

end
#---------------------------------------------------------------------------
# t14_1DStrengthEnvelope/
@testset "t14_1DStrengthEnvelope" begin
    cd(test_dir)
    dir = "t14_1DStrengthEnvelope";
    include(joinpath(dir,"StrengthEnvelop.jl"))

    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=2e-4,atol=1e-10));

    # ---
    # first test runs visco-plastic setup with dt = 10 ka
    @test perform_lamem_test(dir,"1D_VP.dat","1D_VP_Direct_opt",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)
    # ---

    # ---
    # 2nd test runs visco-elasto-plastic setup with dt = 5 ka
    @test perform_lamem_test(dir,"1D_VEP5.dat","1D_VEP5_Direct_opt",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)
    # ---
       
    # ---
    # 3rd test runs visco-elasto-plastic setup with dt = 10 ka
    @test perform_lamem_test(dir,"1D_VEP10.dat","1D_VEP10_Direct_opt",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)
    # ---

    # ---
    # 4th test runs visco-plastic setup with dt = 50 ka
    @test perform_lamem_test(dir,"1D_VEP50.dat","1D_VEP50_Direct_opt",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)
    # ---
    
    # Read output of various simulations:
    VP,  _  = read_LaMEM_timestep("outputVP", 0, dir; last=true);       τII_1 =  Float64.(VP.fields.j2_dev_stress[1,1,:]);
    VEP5,_  = read_LaMEM_timestep("outputVEP5", 0, dir; last=true);     τII_2 =  Float64.(VEP5.fields.j2_dev_stress[1,1,:]);
    VEP10,_ = read_LaMEM_timestep("outputVEP10", 0, dir; last=true);    τII_3 =  Float64.(VEP10.fields.j2_dev_stress[1,1,:]);
    VEP50,_ = read_LaMEM_timestep("outputVEP50", 0, dir; last=true);    τII_4 =  Float64.(VEP50.fields.j2_dev_stress[1,1,:]);

    z       =  VP.z.val[1,1,:]
    phase   =  VP.fields.phase[1,1,:]
    T       =  Float64.(VP.fields.temperature[1,1,:])
    P       =  VP.fields.pressure[1,1,:]*1e6
    τy      =  VP.fields.yield[1,1,:]
    τ_anal  =  Analytical_StrengthEnvelop(phase, T, P, τy)      # analytical solution 

    @test norm(τII_1 - τII_2) ≈ 106.54804495583619
    @test norm(τII_1 - τII_3) ≈ 74.52658531093326
    @test norm(τII_1 - τII_4) ≈ 57.419951900695494
    @test norm(τII_1 - τ_anal) ≈ 148.18742956153582
    
    # Create plot
    if create_plots
        # Plot the strength envelops
        Plot_StrengthEnvelop("StrengthEnvelop_1D.png", dir, z, (τII_1, τII_2, τII_3, τII_4, τ_anal),("Viscoplastic", "VEP dt=5ka", "VEP dt=10ka", "VEP dt=50ka", "Analytical"))
    end

	if clean_files
		clean_test_directory(dir)
	end 

end
#---------------------------------------------------------------------------
@testset "t15_RTI" begin
    dir = "t15_RTI";
    include(joinpath(dir,"RT_analytics.jl"))
    ParamFile = "RTI.dat";
    
    λ       = [0.25, 0.5, 1.0, 1.25, 1.5, 2.0, 4.0]
    q_num   = Compute_RT_growthrate_LaMEM(λ, ParamFile, dir)
    q_anal  = AnalyticalSolution_RTI_FreeSlip(λ)

    @test  norm(q_num - q_anal) ≈ 0.001481104083594891

    # Plot 
    if create_plots
        λ_pl     = range(1e-9,5,100)
        q_anal_pl = AnalyticalSolution_RTI_FreeSlip(λ_pl)
        Plot_growthrate("RTI_analytics_numerics.png", dir, λ,q_num,λ_pl,q_anal_pl)
    end

	if clean_files
		clean_test_directory(dir)
	end 

end
#---------------------------------------------------------------------------
@testset "t16_PhaseTransitions" begin
    cd(test_dir)
    dir = "t16_PhaseTransitions";
    ParamFile = "Plume_PhaseTransitions.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=1e-9));
    
    # Perform tests
    @test perform_lamem_test(dir,ParamFile,"PhaseTransitions",
                            args="-nstep_max 30",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    @test perform_lamem_test(dir,ParamFile,"PhaseTransitions-FreeSlip_p1",
                            args="-open_top_bound 0 -act_press_shift 1 -nstep_max 30", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-9));
    @test perform_lamem_test(dir,"Plume_PhaseTransitions_Melting.dat","PhaseTransitions-Melting_p1",
                            args="-mfmax 0.15",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # Tests phase transitions with X/Z and Box coordinates
    @test perform_lamem_test(dir,"Plume_PhaseTransitions_Box_XZ.dat","PhaseTransitions-XBox",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # Tests phase transition triggered by time       
    @test perform_lamem_test(dir,"TimeTransition.dat","TimeTransition",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
                            
    # Test dike feature using optimized LaMEM
    @test perform_lamem_test(dir,"PhaseTransNotInAirBox_move.dat","PhaseTransNotInAirBox_move",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # Check that it works when one Phase==0; addresses issue #14    
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-9));
    @test perform_lamem_test(dir,"Plume_PhaseTransitions_SwappedPhases.dat","PhaseTransitions-Melting_SwappedPhases_p1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)                 
end
#---------------------------------------------------------------------------
@testset "t17_InflowOutflow" begin
    cd(test_dir)
    dir = "t17_InflowOutflow";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2","|eRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11), (rtol=1e-7,atol=1e-11));
    
    # 2D test
    # InflowOutflow2D_opt
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_2D.dat","InflowOutflow-2D_p1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
    
    # 3D test
    # InflowOutflow3D_opt
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-10), (rtol=1e-4,atol=2e-10));
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_3D.dat","InflowOutflow-3D_p4",
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # Test inflow/outflow conditions in 2D using optimized LaMEM   
    # InflowOutflow2D_Pres_opt 
    acc      = ((rtol=2e-7,atol=2e-7), (rtol=1e-5, atol=1e-6), (rtol=1e-4,atol=2e-8), (rtol=1e-6,atol=1e-9));
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_2D_Perm.dat","InflowOutflow-2D_Perm_p1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # test_3D_Pres():
    # InflowOutflow3D_Pres_opt
    #  keywords = ("|Div|_inf","|Div|_2","|mRes|_2","|eRes|_2")
    acc      = ((rtol=1e-7,atol=1e-7), (rtol=1e-5, atol=1e-5), (rtol=1e-4,atol=1e-8), (rtol=1e-6,atol=1e-9));
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_3D_Perm.dat","InflowOutflow-3D_Perm_p4",
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)         
end
#---------------------------------------------------------------------------
@testset "t18_SimpleShear" begin
    cd(test_dir)
    dir = "t18_SimpleShear";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-8), (rtol=1e-5, atol=2e-8), (rtol=1e-4,atol=2e-5));

    # test_xy
    @test perform_lamem_test(dir,"SimpleShear.dat","SimpleShear_xy",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

end
#---------------------------------------------------------------------------
@testset "t19_CompensatedInflow" begin
    cd(test_dir)
    dir = "t19_CompensatedInflow";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_a
    # t19_CompensatedInflow
    @test perform_lamem_test(dir,"CompensatedInflow_test_2D.dat","CompensatedInflow",
                            args="-nstep_max 10",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # test_b
    # t19_CompensatedInflow3D
    @test perform_lamem_test(dir,"CompensatedInflow_test_3D.dat","CompensatedInflow3D",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # test_migrating ridge
    @test perform_lamem_test(dir,"MigratingRidge_2D.dat","MigratingRidge_2D",
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
    
end
#---------------------------------------------------------------------------
@testset "t20_FSSA" begin
    cd(test_dir)
    dir = "t20_FSSA";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-5,atol=1e-5), (rtol=1e-5, atol=1e-5), (rtol=1e-4,atol=1e-4));

    # t20_FSSA_1_opt
    @test perform_lamem_test(dir,"RTI_FSSA.dat","RTI_FSSA_1",
                            args="-nstep_max 20 -nel_x 50 -nel_z 100",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t21_Passive_Tracer" begin
    cd(test_dir)
    dir = "t21_Passive_Tracer";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-9), (rtol=1e-5, atol=1e-8), (rtol=1e-4,atol=1e-4));

    # test_a
    # t21_Passive_Tracer_Always
    @test perform_lamem_test(dir,"Passive_tracer_ex2D.dat","Passive_tracer-2D_p1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # test_b
    # t21_Passive_Tracer_Condition
    @test perform_lamem_test(dir,"Passive_tracer_ex2D_Condition.dat","Passive_tracer-2D_Condition_p1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t22_RidgeGeom" begin
    cd(test_dir)
    dir = "t22_RidgeGeom";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_2D
    @test perform_lamem_test(dir,"ridge_geom_2D.dat","RidgeGeom2D",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # Test oblique ridge geometry conditions in 3D using optimized LaMEM   
    @test perform_lamem_test(dir,"ridge_geom_oblique_2cores.dat","RidgeGeom_oblique_2cores",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t23_Permeable" begin
    cd(test_dir)
    dir = "t23_Permeable";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=2e-9));

    # test_a
    @test perform_lamem_test(dir,"Permeable.dat","Permeable_p1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t24_Erosion_Sedimentation" begin
    cd(test_dir)
    dir = "t24_Erosion_Sedimentation";
    include(joinpath(dir,"t24_CreateSetup.jl"));      

    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=1e-6), (rtol=1e-5, atol=5e-5), (rtol=2.5e-4,atol=1e-4));
    
    ParamFile = "Erosion_Sedimentation_2D.dat"

    # test_a
    t24_CreateMarkers(dir, ParamFile, NumberCores=2, mpiexec=mpiexec, is64bit=is64bit)
    @test perform_lamem_test(dir,"Erosion_Sedimentation_2D.dat","Erosion_Sedimentation_2D_opt",
                            args="-nstep_max 2",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # test_b
    t24_CreateMarkers(dir, ParamFile, NumberCores=2, mpiexec=mpiexec, is64bit=is64bit)
    @test perform_lamem_test(dir,"Erosion_Sedimentation_2D.dat","Erosion_Sedimentation_2D_deb",
                            args="-nstep_max 2",
                            keywords=keywords, accuracy=acc, cores=2, deb=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t25_APS_Healing" begin
    cd(test_dir)
    dir = "t25_APS_Healing";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-9), (rtol=1e-4,atol=1e-7));

    # test_2D
    @test perform_lamem_test(dir,"APS_Healing2D.dat","APS_Healing2D",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t26_Dike" begin
    cd(test_dir)
    dir = "t26_Dike";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-9), (rtol=1e-5, atol=1e-9), (rtol=1e-4,atol=1e-9));

    # test_M1_2D
    @test perform_lamem_test(dir,"dike_M1_2D.dat","dike_M1_2D",
                            args="-nstep_max 5  -nel_y 2",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # test_M075_2D_2cores
    @test perform_lamem_test(dir,"dike_M075_2D_2cores.dat","dike_M075_2D_2cores",
                            args="-nstep_max 2 -nel_y 2",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # test_variableM
    @test perform_lamem_test(dir,"dike_variableM.dat","dike_variableM",
                            args="-nstep_max 2 -nel_y 2",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # heat_kfac
    keywords = ("|eRes|_2",)
    acc      = ((rtol=1e-4,atol=1e-9),);
    @test perform_lamem_test(dir,"dike_heating_kfac.dat","dike_heating_kfac",
                            args="-nstep_max 2 -nel_y 2",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # heat_rhoA
    keywords = ("|eRes|_2",)
    acc      = ((rtol=1e-4,atol=1e-8),);
    @test perform_lamem_test(dir,"dike_heating_rhoA.dat","dike_heating_rhoA",
                            args="-nstep_max 2 -nel_y 2",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # dyndike_4core.dat
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-9), (rtol=1e-5, atol=1e-9), (rtol=1e-4,atol=1e-9));
    @test perform_lamem_test(dir,"dyndike_4core.dat","dyndike_4core",
                            args="",
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t27_T-dep_Conductivity" begin
    cd(test_dir)
    dir = "t27_T-dep_Conductivity";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_2fields_dike():
    @test perform_lamem_test(dir,"TDep_NuK_Conductivity.dat","TDep_NuK_Conductivity",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t28_HeatRecharge" begin

    cd(test_dir)
    dir = "t28_HeatRecharge";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=5e-7,atol=1e-9), (rtol=1e-6, atol=1e-9), (rtol=2e-5,atol=1e-11));

    # test_recharge1
    @test perform_lamem_test(dir,"FallingBlockHeatReacharge1.dat","HeatRecharge1",
                            args="-nel_x 16 -nel_y 16 -nel_z 16",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # test_recharge2
    acc      = ((rtol=3e-6,atol=5e-6), (rtol=1e-5, atol=1e-5), (rtol=3e-5,atol=2e-5));
    @test perform_lamem_test(dir,"FallingBlockHeatReacharge2.dat","HeatRecharge2",
                            args="-nel_x 16 -nel_y 16 -nel_z 16",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

end
#---------------------------------------------------------------------------
@testset "t29_PermeableSides_VelBoxes" begin
    cd(test_dir)
    dir = "t29_PermeableSides_VelBoxes";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=5e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-5,atol=1e-11));

    # test_permeableSides_VelBoxes
    @test perform_lamem_test(dir,"VelBoxes_Permeable_sides.dat","PermeableSides_VelBoxes",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t30_Timestep_Schedule" begin
    cd(test_dir)
    dir = "t30_Timestep_Schedule";
    
    keywords = ("Actual time step",)
    acc      = ((rtol=1e6,atol=1e-11),);

    # test_TS_Schedule():
    @test perform_lamem_test(dir,"TS_Schedule.dat","TS_Schedule",
                            args="-nel_x 8 -nel_y 8 -nel_z 8",
                            keywords=keywords, accuracy=acc, cores=4, opt=true, split_sign=":", mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t31_geomIO" begin
    cd(test_dir)
    dir = "t31_geomIO";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=2e-3,atol=2e-6), (rtol=5e-3,atol=5e-6), (rtol=5e-3,atol=5e-7));

    # Test if geomIO polygons are read in correctly:
    @test perform_lamem_test(dir,"geomIO_Bulky.dat","geomIO_Bulky",
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)



    # Test if geomIO polygons are read in correctly:
    @test perform_lamem_test(dir,"geomIO_Hollow.dat","geomIO_Hollow",
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

end
#---------------------------------------------------------------------------
@testset "t32_BC_velocity" begin
    cd(test_dir)
    dir = "t32_BC_velocity";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=2e-4,atol=1e-10));

   # Test if boundaries are pushed from front and back inside the model:
    @test perform_lamem_test(dir,"BC_velocity_2D_FB.dat","BC_velocity_2D_FB_opt",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
    
    # Test if boundaries are pushed from left to right and then from right to left:
    @test perform_lamem_test(dir,"BC_velocity_2D_LR.dat","BC_velocity_2D_LR_opt",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t33_Initial_APS" begin
    cd(test_dir)
    dir = "t33_Initial_APS";

    include(joinpath(dir,"initial_aps_setup.jl"))

    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=2e-4,atol=1e-10));
    
    name0 = "no_APS"
    name1 = "APS"
    CreateMarkers_t33(dir, "initial_aps_setup.dat", "./markers_$name0"; NumberCores=1)
    CreateMarkers_t33(dir, "initial_aps_setup.dat", "./markers_$name1"; APS=0.5, NumberCores=1)

    # Test backwards compatibility
    #   Read marker file created without APS column
    #   No passive tracer output
    @test perform_lamem_test(dir,"initial_aps_setup.dat", name0, args="-mark_load_file ./markers_$name0/mdb",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)

    #   APS output on passive tracers works even if no initial APS is set on markers
    @test perform_lamem_test(dir,"initial_aps_setup.dat", name0, args="-mark_load_file ./markers_$name0/mdb -out_ptr 1 -out_ptr_APS 1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)

    # Test initial accumulated plastic strain
    #   Read marker file created with APS=0.5
    #   No passive tracer output
    @test perform_lamem_test(dir,"initial_aps_setup.dat", name1, args="-mark_load_file ./markers_$name1/mdb",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)
    #   APS output on passive tracers
    @test perform_lamem_test(dir,"initial_aps_setup.dat", name1, args="-mark_load_file ./markers_$name1/mdb -out_ptr 1 -out_ptr_APS 1",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=false)

    # Test that plast_strain values are correct in both cases
    #   Includes pvd, marker, and passive tracer output
    mean_APS0, mean_APS1 = compare_APS(dir, "initial_aps_setup.dat", "./markers_$name0", "./markers_$name1")
    #  Verify APS values after 2 timesteps
    @test mean_APS0 == 0.0
    @test mean_APS1 == 0.5

	if clean_files
		clean_test_directory(dir)
	end

end
#---------------------------------------------------------------------------
@testset "t34_3D_2D_push_block" begin
    cd(test_dir)
    dir = "t34_3D_2D_push_block";
    
    ParamFile = "3D_push_block.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=1e-10), (rtol=1e-5,atol=1e-10), (rtol=1e-4,atol=1e-10));
    
    # Test 3D Bezier push block functionality
    @test perform_lamem_test(dir,ParamFile,"3D_push_block_opt", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
							create_expected_file=update_expected, clean_dir=clean_files)

    # Test 2D Bezier push block functionality (push along x-axis only)
    ParamFile = "2D_push_block.dat";
    @test perform_lamem_test(dir,ParamFile,"2D_push_block_opt", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
							create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t35_TopoDiffusion" begin
    cd(test_dir)
    dir = "t35_TopoDiffusion"
    include(joinpath(dir, "TopoDiffusionCreateSetup.jl"))

    keywords = ("|Div|_inf", "|Div|_2", "|mRes|_2")
    acc      = ((rtol=1e-6, atol=1e-6), (rtol=1e-5, atol=5e-5), (rtol=2.5e-4, atol=1e-4))

    ParamFile = "TopoDiffusion.dat"

    TopoDiffusionCreateSetup()

    @test perform_lamem_test(dir, ParamFile, "TopoDiffusion_opt";
        args     = "-nstep_max 3",
        keywords = keywords,
        accuracy = acc,
        cores    = 1,
        opt      = true,
        mpiexec  = mpiexec,
		create_expected_file=update_expected, clean_dir=clean_files)
end

#---------------------------------------------------------------------------
end
#---------------------------------------------------------------------------

