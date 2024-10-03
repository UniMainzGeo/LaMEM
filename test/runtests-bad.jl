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



@testset "t16_PhaseTransitions" begin

    cd(test_dir)
    dir = "t16_PhaseTransitions";
    ParamFile = "Plume_PhaseTransitions.dat";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-6, atol=1e-11), (rtol=2e-6,atol=1e-9));
    
    # Perform tests


    # Test dike feature using optimized LaMEM
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-9));                         
    @test perform_lamem_test(dir,"PhaseTransNotInAirBox_move.dat","PhaseTransNotInAirBox_move.expected",
                            keywords=keywords, accuracy=acc, cores=2, opt=true, mpiexec=mpiexec)
                        
end


@testset "t26_Dike" begin
    cd(test_dir)
    dir = "t26_Dike";
    
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-9), (rtol=1e-5, atol=1e-9), (rtol=1e-4,atol=1e-9));

    # test_variableM
    @test perform_lamem_test(dir,"dike_variableM.dat","dike_variableM.expected",
                            args="-nstep_max 2 -nel_y 2",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)

    # heat_rhoA
    keywords = ("|eRes|_2",)
    acc      = ((rtol=1e-4,atol=1e-8),);
    @test perform_lamem_test(dir,"dike_heating_rhoA.dat","dike_heating_rhoA.expected",
                            args="-nstep_max 2 -nel_y 2",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)

    # dyndike_4core.dat
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-9), (rtol=1e-5, atol=1e-9), (rtol=1e-4,atol=1e-9));
    @test perform_lamem_test(dir,"dyndike_4core.dat","dyndike_4core.expected",
                            args="",
                            keywords=keywords, accuracy=acc, cores=4, opt=true, mpiexec=mpiexec)
    

end


@testset "t27_T-dep_Conductivity" begin
    cd(test_dir)
    dir = "t27_T-dep_Conductivity";
    keywords = ("|Div|_inf","|Div|_2","|mRes|_2")
    acc      = ((rtol=1e-7,atol=1e-11), (rtol=1e-5, atol=1e-11), (rtol=1e-4,atol=1e-11));

    # test_2fields_dike():
    @test perform_lamem_test(dir,"t27_TDep_NuK_Conductivity.dat","t27_TDep_NuK_Conductivity.expected",
                            keywords=keywords, accuracy=acc, cores=1, opt=true, mpiexec=mpiexec)
end



end

