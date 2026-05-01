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






#=


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

