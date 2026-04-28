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

@testset "t10_Compressibility" begin
	


end
#---------------------------------------------------------------------------
end
#---------------------------------------------------------------------------

