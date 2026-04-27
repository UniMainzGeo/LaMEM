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

if "no_superlu" in ARGS
    test_superlu=false
else
    test_superlu=true
end

if "is64bit" in ARGS
    global is64bit=true
else
    global is64bit=false
end

@show use_dynamic_lib test_superlu test_mumps create_plots
include("test_utils.jl")        # test-framework specific functions

test_dir = pwd()

#---------------------------------------------------------------------------

maintenance = true # set to true when designing/debugging tests

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
	# normal mode (do not update anythng and clean all files)
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
    @test perform_lamem_test(dir,ParamFile,"FB1_a_Direct_opt-p1", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    @test perform_lamem_test(dir,ParamFile,"FB1_b_Direct_deb-p1", 
                            keywords=keywords, accuracy=acc, cores=1, deb=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    ParamFile = "FallingBlock_Direct_MUMPS.dat";

    @test perform_lamem_test(dir,ParamFile,"FB1_c_MUMPS_opt-p2", 
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)
end
#---------------------------------------------------------------------------
@testset "t02_FB2_MG" begin
    if test_superlu
        cd(test_dir)
        dir = "t02_FB2_MG";
        
        ParamFile = "FallingBlock_mono_CoupledMG_RedundantCoarse.dat";
        
        keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
        acc      = ((rtol=1e-7,), (rtol=1e-5,), (rtol=1e-4,));
        
        # Perform tests
        @test perform_lamem_test(dir,ParamFile,"FB2_a_CoupledMG_opt-p1", 
                                keywords=keywords, accuracy=acc, cores=4, deb=true, opt=false, mpiexec=mpiexec, debug=false,
                                create_expected_file=update_expected, clean_dir=clean_files)
    end
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

    @test perform_lamem_test(dir,ParamFile,"Sub1_a_Direct_opt-p1", 
                            args="-nstep_max 2",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # Sub1_b_MUMPS_opt                            
    ParamFile = "Subduction_GMG_Particles_MUMPS.dat";

    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=1e-5), (rtol=1e-5,atol=1e-5), (rtol=2.5e-4,atol=1e-3));
        
    CreateMarkers_Subduction(dir, ParamFile, NumberCores=4, mpiexec=mpiexec, is64bit=is64bit)
    @test perform_lamem_test(dir,ParamFile,"Sub1_b_MUMPS_opt-p4", 
                                args="-nstep_max 2",
                                keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                                create_expected_file=update_expected, clean_dir=clean_files)

    # Sub1_c_MUMPS_deb    
    ParamFile = "Subduction_GMG_Particles_MUMPS.dat";
                     
    # writing parallel marker files doesn't work in CI with 64 bit atm 
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-6,atol=2e-6), (rtol=1e-5,atol=3e-6), (rtol=2.5e-4,atol=3e-4));
        
    CreateMarkers_Subduction(dir, ParamFile, NumberCores=4, mpiexec=mpiexec, is64bit=is64bit)
    @test perform_lamem_test(dir,ParamFile,"Sub1_c_MUMPS_deb-p4", 
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
    @test perform_lamem_test(dir,ParamFile,"Sub1_d_MUMPS_MG_VEP_opt-p2", 
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
        @test perform_lamem_test(dir,"localization.dat","Loc1_a_MUMPS_VEP_opt-p4",
                                args="-nstep_max 20", 
                                keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            	create_expected_file=update_expected, clean_dir=clean_files)
    end

    # Loc1_b_MUMPS_VEP_Reg_opt
    @test perform_lamem_test(dir,"localization_eta_min_reg.dat","Loc1_b_MUMPS_VEP_Reg_opt-p4",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)

    # Loc1_c_Direct_VEP_opt                            
    @test perform_lamem_test(dir,"localization.dat","Loc1_c_Direct_VEP_opt-p1",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)


    # Loc1_d_MUMPS_VEP_VPReg_opt
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=5e-2,atol=1e-8), (rtol=2e-3,atol=5e-9), (rtol=2e-3,atol=2e-7));
    
    @test perform_lamem_test(dir,"localization_eta_vp_reg.dat","Loc1_d_MUMPS_VEP_VPReg_opt-p1",
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
    acc      = ((rtol=1e-7,), (rtol=1e-5,), (rtol=1e-2,atol=1e-8));
    
    # Permeability_Direct_opt
    @test perform_lamem_test(dir,ParamFile,"Permeability_direct_opt-p4", 
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
    if test_superlu
        # AdjointGradientScalingLaws_p2
        ParamFile = "RTI_ScalingLaw.dat";
        @test perform_lamem_test(dir,ParamFile,"AdjointGradientScaling_p2",
                                keywords=keywords, accuracy=acc, cores=2, opt=true, split_sign=split_sign, mpiexec=mpiexec,
                                create_expected_file=update_expected, clean_dir=clean_files)
    end
        
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
end
