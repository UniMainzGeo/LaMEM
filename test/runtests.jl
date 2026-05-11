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



	# echo .out file to the screen to see what happens 
	file_content = ""

	open(joinpath(dir,"Erosion_Sedimentation_2D_opt.out"), "r") do io
	    file_content = read(io, String)
	end

	println(file_content)

	

    # test_b
    t24_CreateMarkers(dir, ParamFile, NumberCores=2, mpiexec=mpiexec, is64bit=is64bit)
    @test perform_lamem_test(dir,"Erosion_Sedimentation_2D.dat","Erosion_Sedimentation_2D_deb",
                            args="-nstep_max 2",
                            keywords=keywords, accuracy=acc, cores=2, deb=true, mpiexec=mpiexec,
                            create_expected_file=update_expected, clean_dir=clean_files)


	# echo .out file to the screen to see what happens 
	file_content = ""

	open(joinpath(dir,"Erosion_Sedimentation_2D_deb.out"), "r") do io
	    file_content = read(io, String)
	end

	println(file_content)


end

#---------------------------------------------------------------------------
end
#---------------------------------------------------------------------------

