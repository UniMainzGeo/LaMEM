# This processes all Timestep directories in the current LaMEM 
# directory and creates visualizations of them all

using SpecialFunctions, Printf

# Read LaMEM routines
include("../../scripts/julia/ReadLaMEM_Timestep.jl")

# define a macro that finds all directories with a certain pattern
searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

#using Plots, ColorSchemes
using Plots

FileName    =   "t13_HalfspaceCooling.pvtr";

print("Processing file name $FileName \n")

# Create a new directory that will contain the 2D plots
if ~ispath("Visualization_2D")
    mkdir("Visualization_2D")
end

# Go over all Timestep subdirectories in the current directory
Directories     = searchdir(pwd(),"Timestep_")
iStep           = length(Directories);


DirName     =   Directories[iStep];
print("Processing directory $DirName \n")

# Extract the timestep info from the directory directory name 
id          =   findlast("_",DirName);
Time_Myrs   =   parse(Float64,DirName[id[1]+1:length(DirName)]);

# ------------------------------------------------------------------------------------- 
# Read data from last timestep
data        =   Read_VTR_File(DirName, FileName);
T,x,z       =   ReadField_2D_pVTR(data, "temperature [C]",      "Scalar");     # temperature
# ------------------------------------------------------------------------------------- 

# Compute Analytical solution 
Tstart      = 10;
ThermalAge  = (Tstart+Time_Myrs)*(3600*24*365.25*1e6)
kappa       = 3/(3000*1050);
T_surf      = 20;
T_mantle    = 1350;

T_anal      =   (T_surf-T_mantle).*erfc.((abs.(z).*1e3)./(2*sqrt(kappa*ThermalAge))) .+ T_mantle;

# Create Plot ------------------------------------------------------------------------
plot(T_anal,z,label="Analytics",title="halfspace cooling after $(@sprintf("%.2f", Time_Myrs)) Myrs with 10 Myrs init diffusion")
pl = plot!(T[:,1],z, markershape = :auto,  markersize=1, linewidth = 0, xlabel="Temperature [C]",ylabel="Depth [km]",label="LaMEM")
savefig(pl,"HalfspaceCooling_LaMEM_vs_Analytics.png")