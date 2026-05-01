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
		clean_files = false # no need to clean files, keep working
	end
	
else
	# normal mode (do not update the expected files, just clean the output/work files)
	update_expected = false 
	clean_files     = true
end
#---------------------------------------------------------------------------
@testset "LaMEM Testsuite" verbose=true begin
#---------------------------------------------------------------------------










@testset "t15_RTI" begin
    dir = "t15_RTI";
    include(joinpath(dir,"RT_analytics.jl"))
    ParamFile = "RTI.dat";
    
    λ       = [0.25, 0.5, 1.0, 1.25, 1.5, 2.0, 4.0]
    q_num   = Compute_RT_growthrate_LaMEM(λ, ParamFile, dir)
    q_anal  = AnalyticalSolution_RTI_FreeSlip(λ)

    @test  norm(q_num - q_anal) ≈ 0.001481104 rtol = 1e-4

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

#=
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
	topo_file = "topo.bin"
	
    TopoDiffusionCreateSetup(dir, topo_file)

    @test perform_lamem_test(dir, ParamFile, "TopoDiffusion_opt";
        args     = "-nstep_max 3",
        keywords = keywords,
        accuracy = acc,
        cores    = 1,
        opt      = true,
        mpiexec  = mpiexec,
		create_expected_file=update_expected, clean_dir=clean_files)
	
	if clean_files
		rm(joinpath(dir,topo_file))
	end
end
=#
#---------------------------------------------------------------------------
end
#---------------------------------------------------------------------------

