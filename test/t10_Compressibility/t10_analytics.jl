using LinearAlgebra
using Statistics
using CairoMakie

    
"""
    Sv, Pf, P_hydro, Sh_anal = AnalyticalSolution(rho::Vector{T}, phase_vec::Vector{Integer}, z::Vector{T})

Analytical solution 
"""
function AnalyticalSolution(rho::Vector{T}, phase_vec::Vector{Int64}, z::Vector{T}) where T
    # Define material properties for each of the phases
    # This should obviously be the same as in the input file
    biot = 1.0
    rho_w = 1000  # density of water
    p_fac = [0, 0.1, 0.1, 0.3]  # pore fluid factor
    #rho = rho[:, 1, 1]  # Density
    v_vec = [0.4999, 0.27, 0.4999, 0.27]  # Poissons Ratio

    phase_vec = phase_vec .+1
    nz = length(phase_vec)
    rp = similar(z)
    v = similar(z)

    Sv_anal = similar(z)  # vertical stress
    P_hydro = similar(z)
    Pf_anal = similar(z)  # pore fluid

    # Create arrays with material constants
    for i in reverse(1:nz)
        phase = phase_vec[i]
        v[i] = v_vec[phase_vec[i]]
        rp[i] = p_fac[phase_vec[i]]
        P_hydro[i] = -rho_w * 10 * z[i] * 1e3
    end

    P_hydro[P_hydro .< 0] .= 0

    # Compute vertical stress
    Sv_anal[nz] = 0
    for i in reverse(2:nz)
        rho_mean = (rho[i] + rho[i - 1]) / 2.0
        dz = (z[i] - z[i - 1]) * 1000  # in m
        Sv_anal[i - 1] = Sv_anal[i] + dz * 10 * rho_mean
    end

    Pf_anal = P_hydro + rp .* (Sv_anal - P_hydro)

    Sh_anal = (v ./ (1.0 .- v)) .* (Sv_anal - biot * Pf_anal) + biot * Pf_anal

    return Sv_anal/1e6, Pf_anal/1e6, P_hydro/1e6, Sh_anal/1e6
end


function extract_1D_profiles(data, dir)
    phase_vec = round.(Int64,data.fields.phase[1,1,:]);
    ρ = Float64.(data.fields.density[1,1,:])
    z = data.z.val[1,1,:]

    P = data.fields.pressure
    SHmax_x = data.fields.SHmax[1,:,:,:];
    SHmax_y = data.fields.SHmax[2,:,:,:];
    SHmax_z = data.fields.SHmax[3,:,:,:];

    Txx = data.fields.dev_stress[1,:,:,:];
    Tyy = data.fields.dev_stress[5,:,:,:];
    Tzz = data.fields.dev_stress[9,:,:,:];

    Szz_vec = -(-P[1,1,:] + Tzz[1,1,:]);
    Sxx_vec = -(-P[1,1,:] + Txx[1,1,:]);
    Pf_vec  = data.fields.pore_press[1,1,:];
    τII_vec = data.fields.j2_dev_stress[1,1,:]
    
    return phase_vec,ρ, z, Szz_vec, Sxx_vec, Pf_vec, τII_vec
end

function Plot_vs_analyticalSolution(data, dir, filename="Analytics_vs_LaMEM.png")

    # extract 1D profiles
    phase_vec,ρ, z, Szz_vec, Sxx_vec, Pf_vec, τII_vec = extract_1D_profiles(data, dir)

    # 1D analytical solution
    Sv_a, Pf_a, P_hydro_a, Sh_anal_a = AnalyticalSolution(ρ, phase_vec, z)

    # Open figure 
    f = Figure(resolution = (1500, 800))
    ax = Axis(f[1, 1],  xlabel = "Pressure & Stress [bar]", ylabel = "Depth [km]")
    lines!(ax, Sv_a*10, z,  label = "Analytical σᵥ") 
    lines!(ax, Pf_a*10, z,  label = "Analytical Pf") 
    lines!(ax, P_hydro_a*10, z,  label = "Analytical Ph") 
    lines!(ax, Sh_anal_a*10, z,  label = "Analytical σₕ") 

    scatter!(ax, Pf_vec*10, z,  label = "LaMEM Pf") 
    scatter!(ax, Szz_vec*10, z,  label = "LaMEM σᵥ") 
    scatter!(ax, Sxx_vec*10, z,  label = "LaMEM σₕ") 
    axislegend()

    ax = Axis(f[1, 2],   xlabel = "Stress [MPa]", ylabel = "Depth [km]")
    lines!(ax, τII_vec, z,  label = "LaMEM τII") 
    
    ax = Axis(f[1, 3],  xlabel = "Density [kg/m3]", ylabel = "Depth [km]")
    lines!(ax, ρ, z,  label = "LaMEM τII") 


    save(joinpath(dir,filename), f) 

    return f 
end