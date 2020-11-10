# This plots the reference solution of the original paper along with the solution of the current model
#
# It is handy to see how new code features impact the results

using Plots
using CSV



# Read LaMEM routines
include("../../scripts/julia/ReadLaMEM_Timestep.jl")
include("ReadMinCoord.jl")

# Import reference solution from paper (50)
Analytics = CSV.read("Analytics.csv");
Time_anal = Analytics[:,1];
Zmin_anal = Analytics[:,2];

#--------------------------------------------
# Simulation 1

# run
#run(`rm -rf Timestep*`)
#command_sim = `mpiexec -n 2 ../../bin/opt/LaMEM -ParamFile RTI_FSSA.dat -nstep_max 50 `
#run(command_sim)

Time_vec, Zmin_vec = ReadMinZ("RT_FSSA.pvtr")


plot(Time_anal,Zmin_anal)
plot!(Time_vec, Zmin_vec)
