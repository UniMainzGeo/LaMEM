using SpecialFunctions
import Dates

# Performs LaMEM calculations 
function Compute_RT_growthrate_LaMEM(wav_x, FileName, DirName, OutFile="RTI_test")

    cur_dir = pwd();
    cd(DirName)
    q_num = zero(wav_x)
    for (i,wav) in enumerate(wav_x)
        args = "-coord_x -$(wav/2),$(wav/2) -FreeSurf_Wavelength $wav"
        out = run_lamem_local_test(FileName, 1, args; opt=true, bin_dir="../../bin")  # run LaMEM
        
        data, t = read_LaMEM_timestep(OutFile, 0, pwd(), fields=("velocity [ ]","amplitude [ ]"), surf=true, last=true)   # read surface

        Vz_max = maximum(data.fields.velocity[3])
        A_max  = maximum(data.fields.amplitude)
        q_num[i] = Vz_max/A_max
    end
    cd(cur_dir)

    return q_num
end


"""
    q = AnalyticalSolution_RTI_FreeSlip(λ)

Analytical solution for a RT instability with isoviscous properties, free slip upper/lower bounds and equal layer thicknesses
"""
function AnalyticalSolution_RTI_FreeSlip(λ)
    q_anal = zeros(length(λ))
    for (num,lam) in enumerate(λ)
        k = 2π / lam
        q_anal[num] = ((k^2 + 2) * exp(-k) - exp(-2 * k) - 1.0) / (4 * k * (-2 * k * exp(-k) + exp(-2 * k) - 1))
    end

    return q_anal
end


function Plot_growthrate(filename::String, dir::String, λ,q,λ_anal,q_anal)
    f = Figure(size = (500, 500))
    datetime = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    ax = Axis(f[1, 1],  xlabel = "wavelength [ ]", ylabel = "growthrate [ ]", title=datetime)
    
    lines!(ax, λ_anal,q_anal, label="analyics")
    scatter!(ax, λ,q, label="LaMEM", color=:red)    
    axislegend(position=:rb)

    save(joinpath(dir,filename), f) 

    return nothing
end