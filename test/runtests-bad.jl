# This is the LaMEM testing framework, run through julia
#
# Start it from 
using LaMEM_C
using Test
using GeophysicalModelGenerator
using LaMEM.IO_functions
using CairoMakie
using LaMEM.LaMEM_jll.PETSc_jll

if "use_dynamic_lib" in ARGS
    global use_dynamic_lib=true
else
    global use_dynamic_lib=false
end

test_mumps=true # if we do this later on windows, we have to deactivate this

if "no_superlu" in ARGS
    test_superlu=false
else
    test_superlu=true
end

@show use_dynamic_lib test_superlu test_mumps

test_dir = pwd()

include("test_utils.jl")

# ===================

@testset "LaMEM Testsuite" verbose=true begin



#=

@testset "t4_Localisation" begin
    cd(test_dir)
    dir = "t4_Loc";
    
    ParamFile = "localization.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-10), (rtol=1e-5,atol=2e-9), (rtol=1e-4,atol=1e-7));
    
    # t4_Loc1_d_MUMPS_VEP_VPReg_opt
    @test perform_lamem_test(dir,"localization_eta_vp_reg.dat","t4_Loc1_d_MUMPS_VEP_VPReg_opt-p1.expected",
                            args="-nstep_max 20", 
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)

end


@testset "t12_Temperature_diffusion" begin
    cd(test_dir)
    dir = "t12_Temperature_diffusion";
    include(joinpath(dir,"Temp_setup.jl"))
    ParamFile = "t12_Temperature_diffusion.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=5e-11));

    # ---
    # Perform tests
    CreateMarkers_Temperature(dir, "t12_Temperature_diffusion.dat", "./markers_pT1"; NumberCores=1)

    @test perform_lamem_test(dir,ParamFile,"TpD_a.expected",
                            args="-mark_load_file ./markers_pT1/mdb",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, clean_dir=false,  mpiexec=mpiexec)
    
    data, t1 = Read_LaMEM_timestep("t13", 1, dir); T1=data.fields.temperature[1,1,:]; 
    data, t3 = Read_LaMEM_timestep("t13", 3, dir); T3=data.fields.temperature[1,1,:];
    data, t5 = Read_LaMEM_timestep("t13", 5, dir); T5=data.fields.temperature[1,1,:];
    z = data.z.val[1,1,:]

    T_a5 = Analytical_1D(z, t5)
    @test norm(T_a5 - T5)/length(T5) ≈ 0.03356719876721563

    Plot_Analytics_vs_Numerics(z,T_a5, T5, dir, "T_anal3.png")
    clean_directory(dir)
    # ---
   
    # halfspace cooling test ----
    ParamFile = "t12_Temperature_diffusion_1Dhalfspace.dat"
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2","|T|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=5e-11),  (rtol=1e-1,atol=5e-11));

    @test perform_lamem_test(dir,ParamFile,"t12_Temperature_diffusion-p1.expected",
                args="-printNorms 1",
                keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)

    clean_test_directory(dir)
    # ---
end












@testset "t17_InflowOutflow" begin
    cd(test_dir)
    dir = "t17_InflowOutflow";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2","|eRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11), (rtol=1e-7,atol=1e-11));
    
    # 2D test
    # t17_InflowOutflow2D_opt
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_2D.dat","InflowOutflow-2D_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)
    
    # 3D test
    # t17_InflowOutflow3D_opt
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")


    # Test inflow/outflow conditions in 2D using optimized LaMEM   
    # t17_InflowOutflow2D_Pres_opt 
    acc      = ((rtol=2e-7,atol=2e-7), (rtol=1e-5, atol=1e-6), (rtol=1e-4,atol=2e-8), (rtol=1e-6,atol=1e-9));
    @test perform_lamem_test(dir,"PlumeLithos_Interaction_2D_Perm.dat","InflowOutflow-2D_Perm_p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)

             
end





@testset "t26_Dike" begin
    cd(test_dir)
    dir = "t26_Dike";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-9), (rtol=1e-5, atol=1e-9), (rtol=1e-4,atol=1e-9));


    # test_M075_2D_2cores
    @test perform_lamem_test(dir,"dike_M075_2D_2cores.dat","dike_M075_2D_2cores.expected",
                            args="-nstep_max 2 -nel_y 2",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec)

end









@testset "t14_1DStrengthEnvelope" begin
    cd(test_dir)
    dir = "t14_1DStrengthEnvelope";
    include(joinpath(dir,"StrengthEnvelop.jl"))

    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=2e-4,atol=1e-10));

    # ---
    # first test runs visco-plastic setup with dt = 10 ka
    @test perform_lamem_test(dir,"1D_VP.dat","t14_1D_VP_Direct_opt-p1.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, clean_dir=false, mpiexec=mpiexec)
    # ---

    # ---
    # 2nd test runs visco-elasto-plastic setup with dt = 5 ka
    @test perform_lamem_test(dir,"1D_VEP5.dat","t14_1D_VEP5_Direct_opt-p2.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, clean_dir=false, mpiexec=mpiexec)
    # ---
       
    # ---
    # 3rd test runs visco-elasto-plastic setup with dt = 10 ka
    @test perform_lamem_test(dir,"1D_VEP10.dat","t14_1D_VEP10_Direct_opt-p3.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, clean_dir=false, mpiexec=mpiexec)
    # ---

    # ---
    # 4th test runs visco-plastic setup with dt = 50 ka
    @test perform_lamem_test(dir,"1D_VEP50.dat","t14_1D_VEP50_Direct_opt-p4.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, clean_dir=false, mpiexec=mpiexec)
    # ---
    
    # Read output of various simulations:
    VP,  _  = Read_LaMEM_timestep("outputVP", 0, dir; last=true);       τII_1 =  Float64.(VP.fields.j2_dev_stress[1,1,:]);
    VEP5,_  = Read_LaMEM_timestep("outputVEP5", 0, dir; last=true);     τII_2 =  Float64.(VEP5.fields.j2_dev_stress[1,1,:]);
    VEP10,_ = Read_LaMEM_timestep("outputVEP10", 0, dir; last=true);    τII_3 =  Float64.(VEP10.fields.j2_dev_stress[1,1,:]);
    VEP50,_ = Read_LaMEM_timestep("outputVEP50", 0, dir; last=true);    τII_4 =  Float64.(VEP50.fields.j2_dev_stress[1,1,:]);

    z       =  VP.z.val[1,1,:]
    phase   =  VP.fields.phase[1,1,:]
    T       =  Float64.(VP.fields.temperature[1,1,:])
    P       =  VP.fields.pressure[1,1,:]*1e6
    τy      =  VP.fields.yield[1,1,:]
    τ_anal  =  Analytical_StrengthEnvelop(phase, T, P, τy)      # analytical solution 

    @test norm(τII_1 - τII_2) ≈ 11.4882145f0
    @test norm(τII_1 - τII_3) ≈ 11.532428f0
    @test norm(τII_1 - τII_4) ≈ 13.671384f0
    @test norm(τII_1 - τ_anal) ≈ 147.6532112114033
    
    # Create plot
    Plot_StrengthEnvelop("t14_StrengthEnvelop_1D.png", dir, z, (τII_1, τII_2, τII_3, τII_4, τ_anal),("Viscoplastic", "VEP dt=5ka", "VEP dt=10ka", "VEP dt=50ka", "Analytical"))
    clean_test_directory(dir)
end



=#














end




