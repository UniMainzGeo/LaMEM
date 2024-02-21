# Small routine to compare different plasticity stabilizations in LaMEM, namely:
# 1. No stabilization
# 2. Using the stabilization viscosity per phase (which uses a minimum viscosity per phase)
# 3. Using the viscoplastic regularisation by Duretz et al, which adds a term to the yield stress  (tau_y = tau_y + 2*eta_vp*DIIpl)

# perform runs on the command-line with:
# 1. No stabilization
# ../../bin/opt/LaMEM -ParamFile localization_eta_vp_reg.dat -nel_x 128 -nel_z 64 -eta_vp[1] 1e18 >> no_regularisation.log

# 2. eta_st regularisation (minimum viscosity per phase)
# ../../bin/opt/LaMEM -ParamFile localization_eta_vp_reg.dat -nel_x 128 -nel_z 64 -eta_vp[1] 1e18 -eta_st[1] 1e21 -eta_st[2] 1e21 -eta_st[0] 1e20  >> regularisation_eta_st.log

# 3. eta_vp regularisation (viscoplastic regularisation)
# ../../bin/opt/LaMEM -ParamFile localization_eta_vp_reg.dat -nel_x 128 -nel_z 64 -eta_vp[0] 1e21 -eta_vp[1] 1e21 -eta_vp[2] 1e21 >> regularisation_viscoplastic.log

# 4. eta_min everywhere
# ../../bin/opt/LaMEM -ParamFile localization_eta_vp_reg.dat -nel_x 128 -nel_z 64 -eta_vp[1] 1e18 -eta_st[1] 1e18 -eta_min 1e21 >> regularisation_eta_min.log

# Timings obtained on my machine:
# - no_regularisation.log: 392s
# - regularisation_eta_st.log: 82s
# - regularisation_viscoplastic.log: 57s


#=

# Now lets analyze the logfiles and plot results.
using Makie

function ReadFile(filename, keyword, type=Int64)
    
    f = open(filename)
    lines = readlines(f)
    close(f)
    
    idx = findall(x->occursin(keyword,x), lines)
    values = []
    for i in idx
        push!(values, parse(type, split(lines[i], " ")[end]))
    end

    return @. type(values) 
end

keyword="Number of iterations    :"
iter1 = ReadFile("no_regularisation.log",keyword)
iter2 = ReadFile("regularisation_eta_st.log",keyword)
iter3 = ReadFile("regularisation_viscoplastic.log",keyword)
#iter4 = ReadFile("regularisation_eta_min.log", keyword)

# Create plot
using CairoMakie
x = [1:length(iter1); 1:length(iter2); 1:length(iter3)]
c = [ones(Int64,size(iter1)); ones(Int64,size(iter2))*2; ones(Int64,size(iter4))*3]
iter = [iter1; iter2; iter3]

# Create barplot with results
colors = Makie.wong_colors()
fig = Figure()
ax = Axis(fig[1,1],  title = "localization_eta_vp_reg.dat", xlabel="timestep",ylabel="nonlinear iteration")

barplot!(ax, x, iter, dodge = c, color = colors[c]) # plot
labels = ["no regularisation", "eta_st 1e21", "eta_vp 1e21"]
elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
Legend(fig[1,2], elements, labels) # Legend
fig

save("iter.png", fig, size=(2000,1000))

=#


# Next lets make a few sims with different resolutions and VP regularisation


#../../bin/opt/LaMEM -ParamFile localization_eta_vp_reg.dat -nel_x 64   -nel_z 32 -eta_vp[0] 1e21 -eta_vp[1] 1e21 -eta_vp[2] 1e21 -out_file_name vp_reg_64   
#../../bin/opt/LaMEM -ParamFile localization_eta_vp_reg.dat -nel_x 128 -nel_z  64 -eta_vp[0] 1e21 -eta_vp[1] 1e21 -eta_vp[2] 1e21 -out_file_name vp_reg_128   
#../../bin/opt/LaMEM -ParamFile localization_eta_vp_reg.dat -nel_x 256 -nel_z 128 -eta_vp[0] 1e21 -eta_vp[1] 1e21 -eta_vp[2] 1e21 -out_file_name vp_reg_256   
#../../bin/opt/LaMEM -ParamFile localization_eta_vp_reg.dat -nel_x 512 -nel_z 256 -eta_vp[0] 1e21 -eta_vp[1] 1e21 -eta_vp[2] 1e21 -out_file_name vp_reg_512   

#../../bin/opt/LaMEM -ParamFile localization_eta_vp_reg.dat -nel_x 1024 -nel_z 512 -eta_vp[0] 1e21 -eta_vp[1] 1e21 -eta_vp[2] 1e21 -out_file_name vp_reg_1024 

